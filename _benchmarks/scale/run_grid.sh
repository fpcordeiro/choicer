#!/bin/sh
# The population-scale benchmark grid: two builds of choicer, alternated.
#
#   sh _benchmarks/scale/run_grid.sh <tagA>=<library A> <tagB>=<library B>
#   env: OUT (default _benchmarks/scale/output/grid-<date>), SECTIONS (default
#   "A B C D E F G H I J"), MC_SHIM (the malloc shim, macOS: memory runs through
#   output/Rmc, see make_rmc.sh; without it, time runs only)
#
# Each library holds one choicer build (R CMD INSTALL --library=<dir>); every
# run asserts it loaded choicer from there. Short configurations: three rounds
# of time runs, alternated AB, BA, AB, then one memory run per build under the
# shim. Long ones (10^8 rows, minutes per call): one time run per build. Each
# run waits (up to 10 minutes) while the machine's processes use more than 3
# cores, and the log records the load, CPU, available memory, swap and top
# processes; 10^8-row runs need enough available memory, and a run during which
# swap grew by more than 256 MB is renamed *.swapped.rds. Builds older than
# 0.3.0 lack some kernels (the clustered meat, the layout probe), so those
# configurations run on B only. Results: one .rds per run (bench_kernels.R),
# prep_grid.csv and prep_time.csv (prep_mem.R), run.log.
#
# Sections: A 10^8-row designs (gradient at 1 and 11 threads, tastes,
# prediction, BHHH and meats on B, store-mode prediction, the store-mode cube);
# B the gradient grid and store mode; C heavy panels (S = 500, 1000); D the
# Hessian, BHHH, scores, meats, tastes and prediction; E data preparation;
# F the per-call layout (B); G B's Hessian at 10^8 rows; H the multinomial and
# nested logits' derivative kernels and meats (bench_mnl_nl.R); I their fits
# per se_method (fit_mnl_nl.R, one round each); J the multinomial and nested
# logits' per-situation kernels (gradient, BHHH, prediction, shares, BLP),
# every model's elasticity and diversion kernels, their heap, and instruction
# counts per situation (bench_pmc.R, macOS). In H, kernels a build lacks
# (the clustered meat before 0.3.0) record their error, and the dense kernels
# of builds before 0.3.0 run once (at 17,000 parameters on a few hundred
# situations: a full call would take days).
D=$(cd "$(dirname "$0")" && pwd)
A=${1%%=*}; LIB_A=${1#*=}; B=${2%%=*}; LIB_B=${2#*=}
[ -n "$A" ] && [ -n "$B" ] && [ -d "$LIB_A" ] && [ -d "$LIB_B" ] || {
  echo "usage: sh run_grid.sh <tagA>=<libA> <tagB>=<libB>"; exit 1; }
OUT=${OUT:-$D/output/grid-$(date +%Y%m%d-%H%M%S)}; mkdir -p "$OUT"
LOG=$OUT/run.log
SECTIONS=${SECTIONS:-"A B C D E F G H I J"}
SCRIPT=
ROUND=

lib_of() { [ "$1" = "$A" ] && echo "$LIB_A" || echo "$LIB_B"; }
ps_top() { ps -A -o %cpu=,comm= | sort -rn | head -4 | awk '{n=split($2,p,"/"); printf "%s%% %s; ", $1, p[n]}'; }
cpu_all() { ps -A -o %cpu= | awk '{s+=$1} END {printf "%d", s}'; }
load_avg() { if [ -r /proc/loadavg ]; then cut -d' ' -f1-3 /proc/loadavg; else sysctl -n vm.loadavg; fi; }
avail_gb() {
  if [ -r /proc/meminfo ]; then awk '/^MemAvailable:/ {printf "%.1f", $2 * 1024 / 1e9}' /proc/meminfo
  else vm_stat | awk '/page size of/ {ps=$8} /Pages free/ {f=$3} /Pages inactive/ {i=$3}
    /Pages speculative/ {s=$3} /Pages purgeable/ {p=$3}
    END {gsub(/\./,"",f); gsub(/\./,"",i); gsub(/\./,"",s); gsub(/\./,"",p); printf "%.1f", (f+i+s+p)*ps/1e9}'; fi; }
swap_mb() {
  if [ -r /proc/meminfo ]; then awk '/^SwapTotal:/ {t=$2} /^SwapFree:/ {f=$2} END {printf "%d", (t-f)/1024}' /proc/meminfo
  else sysctl -n vm.swapusage | awk '{v=$6; gsub(/M/,"",v); printf "%d", v}'; fi; }
wait_quiet() {
  n=0
  while [ "$(cpu_all)" -gt 300 ] && [ $n -lt 60 ]; do sleep 10; n=$((n+1)); done
  [ $n -ge 60 ] && echo "  BUSY: still $(cpu_all)% after 10 min ($(ps_top))" >> "$LOG"
  return 0
}
mem_ok() {  # <GB needed> <label>
  av=$(avail_gb)
  awk -v a="$av" -v n="$1" 'BEGIN {exit !(a >= n)}' && return 0
  echo "SKIP $2: available ${av} GB < $1 GB" >> "$LOG"; return 1
}
run1() {  # <tag> <time|mem> <name> args...
  tag=$1; mode=$2; name=$3; shift 3
  f=$OUT/${name}_${tag}_${mode}${ROUND:+_r$ROUND}.rds
  [ -f "$f" ] && return 0
  if [ "$mode" = mem ] && [ -z "$MC_SHIM" ]; then return 0; fi
  [ "$mode" = time ] && wait_quiet
  s0=$(swap_mb)
  echo "---- $(date +%H:%M:%S) $tag $mode $name r${ROUND:-1} | $* | load $(load_avg) | cpu $(cpu_all)% | avail $(avail_gb) GB | swap ${s0}M | $(ps_top)" >> "$LOG"
  script=${SCRIPT:-bench_kernels.R}
  if [ "$mode" = mem ]; then
    env R_LIBS="$(lib_of $tag)" MC_SHIM="$MC_SHIM" CHOICER_SCALE_OUT="$OUT" "$D/output/Rmc" --no-echo --no-restore \
      --file="$D/$script" --args tag=$tag mode=mem out="$f" "$@" >> "$LOG" 2>&1
  else
    env R_LIBS="$(lib_of $tag)" CHOICER_SCALE_OUT="$OUT" Rscript --no-echo --no-restore "$D/$script" \
      tag=$tag mode=time out="$f" "$@" >> "$LOG" 2>&1
  fi
  st=$?; s1=$(swap_mb)
  echo "     end $(date +%H:%M:%S) exit $st | swap ${s1}M" >> "$LOG"
  if [ $((s1 - s0)) -gt 256 ] && [ -f "$f" ]; then mv "$f" "${f%.rds}.swapped.rds"; echo "  SWAPPED: +$((s1 - s0)) MB" >> "$LOG"; fi
  return 0
}
pair3() { name=$1; shift
  for r in 1 2 3; do ROUND=$r
    if [ $r = 2 ]; then run1 $B time $name "$@"; run1 $A time $name "$@"
    else run1 $A time $name "$@"; run1 $B time $name "$@"; fi
  done; ROUND=; }
pair6() { name=$1; shift
  for r in 1 2 3 4 5 6; do ROUND=$r
    if [ $((r % 2)) = 0 ]; then run1 $B time $name "$@"; run1 $A time $name "$@"
    else run1 $A time $name "$@"; run1 $B time $name "$@"; fi
  done; ROUND=; }
once() { name=$1; shift; ROUND=1; run1 $A time $name "$@"; run1 $B time $name "$@"; ROUND=; }
once_b() { name=$1; shift; ROUND=1; run1 $B time $name "$@"; ROUND=; }
mem2() { name=$1; shift; run1 $A mem $name "$@"; run1 $B mem $name "$@"; }
has() { case " $SECTIONS " in *" $1 "*) return 0;; esac; return 1; }

echo "==== start $(date) | $A=$LIB_A $B=$LIB_B | load $(load_avg) | avail $(avail_gb) GB | swap $(swap_mb)M" >> "$LOG"
C2M="regime=census U=2000000 S=100"
CL150="regime=claims U=150000 S=100"
if has A; then
  mem_ok 14 census2m && once census2m_grad $C2M threads=1,11 reps=2 what=grad
  mem_ok 13 claims150k && once claims150k_grad $CL150 threads=1,11 reps=2 what=grad
  once small2500k_grad regime=small U=2500000 S=100 threads=1,11 reps=2 what=grad
  mem_ok 14 census2m && once census2m_tastes $C2M threads=11 reps=1 what=tastes
  mem_ok 14 census2m && once_b census2m_main $C2M threads=11 reps=1 what=predict,bhhh,meat_cluster
  mem_ok 26 census2m-A-predict && { ROUND=1; run1 $A time census2m_predict $C2M threads=11 reps=1 what=predict; ROUND=; }
  mem_ok 14 census2m && once_b census2m_storepred $C2M threads=11 reps=1 draws=store what=predict
  mem_ok 21 census2m-store && once_b census2m_store $C2M threads=11 reps=1 draws=store what=grad
  mem_ok 14 census2m && mem2 census2m_grad $C2M threads=11 what=grad,tastes
  mem_ok 14 census2m && run1 $B mem census2m_main $C2M threads=11 what=predict,bhhh,meat_cluster
  mem_ok 13 claims150k && once claims150k_other $CL150 threads=11 reps=1 what=bhhh,scores,tastes
  mem_ok 13 claims150k && once_b claims150k_predict $CL150 threads=11 reps=1 what=predict
  mem_ok 18 claims150k-A-predict && { ROUND=1; run1 $A time claims150k_predict $CL150 threads=11 reps=1 what=predict; ROUND=; }
  mem_ok 13 claims150k && mem2 claims150k_post $CL150 threads=11 what=grad,bhhh,scores,tastes
  mem_ok 13 claims150k && run1 $B mem claims150k_predict $CL150 threads=11 what=predict
  mem_ok 18 claims150k-A-predict && run1 $A mem claims150k_predict $CL150 threads=11 what=predict
  mem2 small2500k_grad regime=small U=2500000 S=100 threads=11 what=grad
fi
if has B; then
  pair3 claims1500_grad regime=claims U=1500 S=100 threads=1,2,5,11 reps=5 what=grad
  pair3 census20k_grad regime=census U=20000 S=100 threads=1,2,5,11 reps=5 what=grad
  pair3 claims15k_grad regime=claims U=15000 S=100 threads=1,2,5,11 reps=3 what=grad
  pair3 census200k_grad regime=census U=200000 S=100 threads=1,2,5,11 reps=3 what=grad
  pair3 small250k_grad regime=small U=250000 S=100 threads=1,11 reps=3 what=grad
  for c in claims:1500 census:20000 claims:15000 census:200000 small:250000; do
    mem2 ${c%%:*}${c##*:}_grad regime=${c%%:*} U=${c##*:} S=100 threads=11 what=grad
  done
  pair3 small250k_store regime=small U=250000 S=100 threads=1,11 reps=3 draws=store what=grad
  pair3 claims15k_store regime=claims U=15000 S=100 threads=11 reps=3 draws=store what=grad,predict
  pair3 census200k_store regime=census U=200000 S=100 threads=11 reps=3 draws=store what=grad,predict
  mem2 claims15k_store regime=claims U=15000 S=100 threads=11 draws=store what=grad,predict
  mem2 census200k_store regime=census U=200000 S=100 threads=11 draws=store what=grad,predict
fi
if has C; then
  pair3 claims1500_S1000 regime=claims U=1500 S=1000 threads=1,11 reps=2 what=grad,tastes
  pair3 claims15k_S500 regime=claims U=15000 S=500 threads=11 reps=2 what=grad
  pair3 claims15k_S1000 regime=claims U=15000 S=1000 threads=11 reps=1 what=grad
  pair3 census200k_S1000 regime=census U=200000 S=1000 threads=11 reps=1 what=grad
  mem2 claims1500_S1000 regime=claims U=1500 S=1000 threads=11 what=grad,tastes
  mem2 claims15k_S1000 regime=claims U=15000 S=1000 threads=11 what=grad
  mem2 census200k_S1000 regime=census U=200000 S=1000 threads=11 what=grad
fi
if has D; then
  pair3 claims1500_other regime=claims U=1500 S=100 threads=11 reps=5 what=bhhh,scores,tastes,predict
  once claims1500_hess regime=claims U=1500 S=100 threads=1,11 reps=1 warm=0 what=hess
  pair3 census2k_hess regime=census U=2000 S=100 threads=11 reps=1 warm=0 what=hess
  pair3 census20k_other regime=census U=20000 S=100 threads=11 reps=3 what=scores,tastes,predict
  pair3 census20k_bhhh regime=census U=20000 S=100 threads=11 reps=1 warm=0 what=bhhh
  once census20k_hess regime=census U=20000 S=100 threads=11 reps=1 warm=0 what=hess
  pair3 claims15k_other regime=claims U=15000 S=100 threads=11 reps=3 what=bhhh,scores,tastes,predict
  once claims15k_hess regime=claims U=15000 S=100 threads=11 reps=1 warm=0 what=hess
  pair3 census200k_other regime=census U=200000 S=100 threads=11 reps=3 what=tastes,predict
  for r in 1 2 3; do ROUND=$r; run1 $B time census200k_main regime=census U=200000 S=100 threads=11 reps=3 what=bhhh,meat_cluster; done; ROUND=
  ROUND=1; run1 $A time census200k_meatr regime=census U=200000 S=100 threads=11 reps=1 warm=0 what=meatr_robust,meatr_cluster; ROUND=
  once_b census200k_hess regime=census U=200000 S=100 threads=11 reps=1 warm=0 what=hess
  mem2 claims1500_other regime=claims U=1500 S=100 threads=11 what=hess,bhhh,scores,tastes,predict
  mem2 census2k_hess regime=census U=2000 S=100 threads=11 what=hess
  mem2 census20k_other regime=census U=20000 S=100 threads=11 what=bhhh,scores,tastes,predict
  mem2 claims15k_other regime=claims U=15000 S=100 threads=11 what=bhhh,scores,tastes,predict
  mem2 census200k_other regime=census U=200000 S=100 threads=11 what=tastes,predict
  run1 $B mem census200k_main regime=census U=200000 S=100 threads=11 what=bhhh,meat_cluster,hess
  run1 $A mem census200k_meatr regime=census U=200000 S=100 threads=11 what=meatr_robust,meatr_cluster
fi
if has E; then
  echo "---- $(date +%H:%M:%S) data preparation | load $(load_avg) | $(ps_top)" >> "$LOG"
  for round in 1 2 3; do for spec in claims:15000 census:200000; do
    if [ $round = 2 ]; then order="$B $A"; else order="$A $B"; fi
    for model in mxl mnl nl; do for t in $order; do
      [ -n "$MC_SHIM" ] && env R_LIBS="$(lib_of $t)" MC_SHIM="$MC_SHIM" CHOICER_SCALE_OUT="$OUT" "$D/output/Rmc" \
        --no-echo --no-restore --file="$D/prep_mem.R" --args regime=${spec%%:*} U=${spec##*:} model=$model \
        tag=$t csv="$OUT/prep_grid.csv" >> "$LOG" 2>&1
    done; done
    for t in $order; do wait_quiet
      env R_LIBS="$(lib_of $t)" CHOICER_SCALE_OUT="$OUT" Rscript --no-echo --no-restore "$D/prep_mem.R" \
        regime=${spec%%:*} U=${spec##*:} model=mxl tag=$t csv="$OUT/prep_time.csv" >> "$LOG" 2>&1
    done
  done; done
fi
if has F; then
  once_b small2500k_layout regime=small U=2500000 S=100 Kx=1 Kw=1 threads=1,11 reps=5 what=layout,grad
  once_b small2500k_layout25 regime=small U=2500000 S=25 Kx=1 Kw=1 threads=11 reps=3 what=grad
  mem_ok 6 census2m-K1 && once_b census2m_layout regime=census U=2000000 S=100 Kx=1 Kw=1 threads=1,11 reps=5 what=layout,grad
  mem_ok 6 census2m-K1 && once_b census2m_layout25 regime=census U=2000000 S=25 Kx=1 Kw=1 threads=11 reps=3 what=grad
fi
if has G; then
  mem_ok 14 census2m && once_b census2m_hess $C2M threads=11 reps=1 warm=0 what=hess
fi
if has H; then
  SCRIPT=bench_mnl_nl.R
  for m in mnl nl; do
    # Variances on census-style and claims-style designs; main's dense kernels
    # once per configuration where a call takes minutes.
    pair3 ${m}_census20k_var model=$m regime=census U=20000 threads=1,5 reps=1 warm=0 what=hess,bhhh,robust
    pair3 ${m}_census20k_meat model=$m regime=census U=20000 threads=1,5 reps=1 warm=0 what=meat_cluster,meatr_robust,meatr_cluster
    pair3 ${m}_census20k_meatG5 model=$m regime=census U=20000 clusters=5 threads=5 reps=1 warm=0 what=meat_cluster,meatr_cluster
    pair3 ${m}_claims1500_var model=$m regime=claims U=1500 threads=1,5 reps=3 what=hess,bhhh,robust,meat_cluster,meatr_cluster
    pair3 ${m}_claims15k_var model=$m regime=claims U=15000 threads=1,5 reps=1 warm=0 what=hess,bhhh,robust,meat_cluster,meatr_cluster
    once ${m}_census200k_bhhh model=$m regime=census U=200000 threads=11 reps=1 warm=0 what=bhhh
    for r in 1 2 3; do ROUND=$r; run1 $B time ${m}_census200k_main model=$m regime=census U=200000 threads=1,11 reps=3 what=hess,bhhh,robust,meat_cluster; done; ROUND=
    once ${m}_census20k17k_sub model=$m regime=census20k U=20000 situations=500 threads=5 reps=1 warm=0 what=hess,bhhh
    once_b ${m}_census20k17k model=$m regime=census20k U=20000 threads=5,11 reps=1 warm=0 what=hess,bhhh,robust,meat_cluster
    # No regression where the compact block is the dense one (4 alternatives,
    # for the nested logit in two nests of two), and the gradient, unchanged
    # by Phase 8A.
    pair6 ${m}_small2500k model=$m regime=small U=2500000 nests=2 threads=11 reps=3 what=grad,hess,bhhh
    pair3 ${m}_small2500k_1t model=$m regime=small U=2500000 nests=2 threads=1 reps=3 what=grad,hess,bhhh
    mem2 ${m}_census20k_var model=$m regime=census U=20000 threads=5 what=hess,bhhh,robust,meatr_robust,meatr_cluster
    run1 $B mem ${m}_census200k_main model=$m regime=census U=200000 threads=11 what=hess,bhhh,robust,meat_cluster
  done
  SCRIPT=
fi
if has I; then
  for spec in census:20000 claims:15000; do
    for m in mnl nl; do for t in $A $B; do
      f=$OUT/${m}_fit_${spec%%:*}${spec##*:}_${t}.rds
      [ -f "$f" ] && continue
      wait_quiet
      echo "---- $(date +%H:%M:%S) $t fit $m ${spec} | load $(load_avg) | cpu $(cpu_all)% | $(ps_top)" >> "$LOG"
      env R_LIBS="$(lib_of $t)" CHOICER_SCALE_OUT="$OUT" Rscript --no-echo --no-restore "$D/fit_mnl_nl.R" \
        tag=$t model=$m regime=${spec%%:*} U=${spec##*:} threads=5 out="$f" >> "$LOG" 2>&1
    done; done
  done
fi
if has J; then
  SCRIPT=bench_mnl_nl.R
  for m in mnl nl; do
    # Per-situation work where it shows (4 alternatives; NL two nests of two)
    pair3 ${m}_small2500k_1t_k model=$m regime=small U=2500000 nests=2 threads=1 reps=3 \
      what=grad,grad_opg,bhhh,hess,predict,shares,elas,dr,blp
    pair3 ${m}_small2500k_5t_k model=$m regime=small U=2500000 nests=2 threads=5 reps=3 \
      what=grad,grad_opg,bhhh,hess,predict,shares,elas,dr,blp
    pair6 ${m}_small2500k_11t_k model=$m regime=small U=2500000 nests=2 threads=11 reps=3 \
      what=grad,grad_opg,bhhh,hess,predict,shares,elas,dr,blp
    # Choice sets of 10-100 and 20-200
    pair3 ${m}_census20k_k model=$m regime=census U=20000 threads=1,5 reps=3 \
      what=grad,grad_opg,bhhh,hess,predict,shares,elas,dr,blp
    pair3 ${m}_claims15k_k model=$m regime=claims U=15000 threads=1,5 reps=3 \
      what=grad,grad_opg,bhhh,hess,predict,shares,elas,dr,blp
    # J = 20,000 elasticities and diversion: per-thread copies (B's budget)
    # against one shared result (acc=0), with the outside option
    once_b ${m}_census20kJ_jj model=$m regime=census20k U=20000 threads=5,11 reps=1 warm=0 ioo=1 what=elas,dr
    once_b ${m}_census20kJ_jj_shared model=$m regime=census20k U=20000 threads=5,11 reps=1 warm=0 ioo=1 acc=0 what=elas,dr
    mem2 ${m}_small250k_k model=$m regime=small U=250000 nests=2 threads=1,11 \
      what=grad,grad_opg,bhhh,hess,scores,meat_cluster,predict,shares,elas,dr,blp
    mem2 ${m}_census20k_k model=$m regime=census U=20000 threads=1,11 \
      what=grad,grad_opg,bhhh,hess,scores,meat_cluster,predict,shares,elas,dr,blp
  done
  SCRIPT=
  pair3 mxl_census20k_jj regime=census U=20000 S=100 threads=1,5 reps=1 what=elas,elas_r,dr,dr_r
  mem2 mxl_census20k_jj regime=census U=20000 S=100 threads=1,11 what=elas,dr
  # Instructions and cycles per situation at one thread, both builds alternated
  # (the mixed logit on 100,000 decision makers: its store-mode cube)
  for r in 1 2 3; do for t in $A $B; do
    for spec in mnl,nl:2500000 mxl:100000; do
      f=$OUT/pmc_small_${spec%%:*}_${t}_r$r.rds
      [ -f "$f" ] && continue
      wait_quiet
      echo "---- $(date +%H:%M:%S) $t pmc ${spec} r$r | load $(load_avg) | $(ps_top)" >> "$LOG"
      env KMP_BLOCKTIME=0 R_LIBS="$(lib_of $t)" CHOICER_SCALE_OUT="$OUT" Rscript --no-echo --no-restore \
        "$D/bench_pmc.R" tag=$t models=${spec%%:*} kernels=grad,grad_opg,bhhh,hess U=${spec##*:} S=100 \
        out="$f" >> "$LOG" 2>&1
    done
  done; done
fi
echo "==== done $(date) | load $(load_avg)" >> "$LOG"
