// Allocation counter + RSS sampler for measuring choicer kernels (macOS).
//
// Load with DYLD_INSERT_LIBRARIES=<path>/memcount.dylib so malloc & co. are
// interposed in every image of the R process; then dyn.load() the same file
// from R and use .C():
//   .C("mc_reset")                       zero the counters, start counting
//   .C("mc_get", out = double(6))        n_alloc, bytes_alloc, live_now,
//                                        peak_live (since reset), n_free,
//                                        bytes_freed
//   .C("mc_rss_start", as.double(ms))    start a sampler thread (period ms)
//   .C("mc_rss_stop", out = double(4))   stop it: peak resident, peak
//                                        footprint, resident now, footprint now
// Counting is off until the first mc_reset(). Live bytes are net of frees of
// blocks allocated before the reset, so peak_live is the high-water mark of
// net heap growth since the reset.
#include <malloc/malloc.h>
#include <mach/mach.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct {
  const void *replacement;
  const void *replacee;
} interpose_t;
#define DYLD_INTERPOSE(_replacement, _replacee)                                \
  __attribute__((used)) static const interpose_t _interpose_##_replacee        \
      __attribute__((section("__DATA,__interpose"))) = {                       \
          (const void *)(unsigned long)&_replacement,                          \
          (const void *)(unsigned long)&_replacee};

static _Atomic int64_t n_alloc = 0, bytes_alloc = 0, live = 0, peak = 0;
static _Atomic int64_t n_free = 0, bytes_freed = 0;
static _Atomic int enabled = 0;

static inline void on_alloc(void *p) {
  if (!p || !atomic_load_explicit(&enabled, memory_order_relaxed)) return;
  const int64_t sz = (int64_t)malloc_size(p);
  atomic_fetch_add_explicit(&n_alloc, 1, memory_order_relaxed);
  atomic_fetch_add_explicit(&bytes_alloc, sz, memory_order_relaxed);
  const int64_t l =
      atomic_fetch_add_explicit(&live, sz, memory_order_relaxed) + sz;
  int64_t pk = atomic_load_explicit(&peak, memory_order_relaxed);
  while (l > pk &&
         !atomic_compare_exchange_weak_explicit(
             &peak, &pk, l, memory_order_relaxed, memory_order_relaxed)) {
  }
}

static inline void on_free_size(int64_t sz) {
  if (!atomic_load_explicit(&enabled, memory_order_relaxed)) return;
  atomic_fetch_add_explicit(&n_free, 1, memory_order_relaxed);
  atomic_fetch_add_explicit(&bytes_freed, sz, memory_order_relaxed);
  atomic_fetch_sub_explicit(&live, sz, memory_order_relaxed);
}

static void *mc_malloc(size_t s) {
  void *p = malloc(s);
  on_alloc(p);
  return p;
}
static void *mc_calloc(size_t n, size_t s) {
  void *p = calloc(n, s);
  on_alloc(p);
  return p;
}
static void *mc_valloc(size_t s) {
  void *p = valloc(s);
  on_alloc(p);
  return p;
}
static void *mc_realloc(void *old, size_t s) {
  const int64_t old_sz = old ? (int64_t)malloc_size(old) : 0;
  void *p = realloc(old, s);
  if (p) {
    if (old) on_free_size(old_sz);
    on_alloc(p);
  }
  return p;
}
static void mc_free(void *p) {
  if (p) on_free_size((int64_t)malloc_size(p));
  free(p);
}
static int mc_posix_memalign(void **pp, size_t a, size_t s) {
  const int r = posix_memalign(pp, a, s);
  if (r == 0) on_alloc(*pp);
  return r;
}
static void *mc_aligned_alloc(size_t a, size_t s) {
  void *p = aligned_alloc(a, s);
  on_alloc(p);
  return p;
}

DYLD_INTERPOSE(mc_malloc, malloc)
DYLD_INTERPOSE(mc_calloc, calloc)
DYLD_INTERPOSE(mc_valloc, valloc)
DYLD_INTERPOSE(mc_realloc, realloc)
DYLD_INTERPOSE(mc_free, free)
DYLD_INTERPOSE(mc_posix_memalign, posix_memalign)
DYLD_INTERPOSE(mc_aligned_alloc, aligned_alloc)

void mc_reset(void) {
  atomic_store(&enabled, 0);
  atomic_store(&n_alloc, 0);
  atomic_store(&bytes_alloc, 0);
  atomic_store(&live, 0);
  atomic_store(&peak, 0);
  atomic_store(&n_free, 0);
  atomic_store(&bytes_freed, 0);
  atomic_store(&enabled, 1);
}

void mc_get(double *out) {
  out[0] = (double)atomic_load(&n_alloc);
  out[1] = (double)atomic_load(&bytes_alloc);
  out[2] = (double)atomic_load(&live);
  out[3] = (double)atomic_load(&peak);
  out[4] = (double)atomic_load(&n_free);
  out[5] = (double)atomic_load(&bytes_freed);
}

// --- RSS / footprint sampler -------------------------------------------------
static pthread_t sampler;
static _Atomic int sampling = 0;
static _Atomic int64_t peak_rss = 0, peak_fp = 0;
static useconds_t period_us = 1000;

static void read_mem(int64_t *rss, int64_t *fp) {
  task_vm_info_data_t info;
  mach_msg_type_number_t count = TASK_VM_INFO_COUNT;
  if (task_info(mach_task_self(), TASK_VM_INFO, (task_info_t)&info, &count) ==
      KERN_SUCCESS) {
    *rss = (int64_t)info.resident_size;
    *fp = (int64_t)info.phys_footprint;
  } else {
    *rss = *fp = -1;
  }
}

static void *sampler_main(void *arg) {
  (void)arg;
  while (atomic_load(&sampling)) {
    int64_t rss, fp;
    read_mem(&rss, &fp);
    if (rss > atomic_load(&peak_rss)) atomic_store(&peak_rss, rss);
    if (fp > atomic_load(&peak_fp)) atomic_store(&peak_fp, fp);
    usleep(period_us);
  }
  return NULL;
}

void mc_rss_start(double *ms) {
  if (atomic_load(&sampling)) return;
  period_us = (useconds_t)(ms[0] * 1000.0);
  if (period_us < 100) period_us = 100;
  int64_t rss, fp;
  read_mem(&rss, &fp);
  atomic_store(&peak_rss, rss);
  atomic_store(&peak_fp, fp);
  atomic_store(&sampling, 1);
  pthread_create(&sampler, NULL, sampler_main, NULL);
}

void mc_rss_stop(double *out) {
  if (atomic_load(&sampling)) {
    atomic_store(&sampling, 0);
    pthread_join(sampler, NULL);
  }
  int64_t rss, fp;
  read_mem(&rss, &fp);
  out[0] = (double)atomic_load(&peak_rss);
  out[1] = (double)atomic_load(&peak_fp);
  out[2] = (double)rss;
  out[3] = (double)fp;
}

void mc_rss_now(double *out) {
  int64_t rss, fp;
  read_mem(&rss, &fp);
  out[0] = (double)rss;
  out[1] = (double)fp;
}

void mc_disable(void) { atomic_store(&enabled, 0); }
