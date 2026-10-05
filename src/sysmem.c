/* Physical memory of the machine, for the mixed logit kernels' memory check
   (src/mxlogit.cpp). In C and without R's headers, which clash with
   <windows.h>. */
#if defined(__APPLE__) && !defined(_DARWIN_C_SOURCE)
#define _DARWIN_C_SOURCE /* _SC_PHYS_PAGES, hidden under strict POSIX macros */
#endif
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

/* Bytes of physical memory, or 0 when they cannot be read. */
double choicer_physical_memory(void);
double choicer_physical_memory(void) {
#ifdef _WIN32
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  return GlobalMemoryStatusEx(&status) ? (double) status.ullTotalPhys : 0.0;
#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
  const long pages = sysconf(_SC_PHYS_PAGES);
  const long page_size = sysconf(_SC_PAGESIZE);
  return pages > 0 && page_size > 0 ? (double) pages * (double) page_size : 0.0;
#else
  return 0.0;
#endif
}
