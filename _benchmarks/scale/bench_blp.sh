#!/bin/sh
# The share inversion's benchmark (bench_blp.R) on two builds, alternated.
#
#   sh _benchmarks/scale/bench_blp.sh <tagA>=<library A> <tagB>=<library B> [<configs>]
#   env: OUT (default _benchmarks/scale/output/blp-<date>), ROUNDS (3), QUIET_CPU
#   (300: the CPU percentage of all processes above which a run waits)
#
# <configs> is a file with one configuration per line, bench_blp.R's options
# (lines starting with # are skipped); the default is the grid below: every
# method from the fitted constants on claims, census and small for each
# model, at 1 and 11 threads. Each configuration runs ROUNDS rounds, the
# builds alternated (AB, BA, AB, ...), one R process per run, after waiting
# (up to 10 minutes) while the machine's processes use more than QUIET_CPU
# percent of a core (300 by default); the
# log records the load and the top processes of every run. Results: one .rds
# per run, named <tag>_<line>_<round>.rds, and run.log.
D=$(cd "$(dirname "$0")" && pwd)
A=${1%%=*}; LIB_A=${1#*=}; B=${2%%=*}; LIB_B=${2#*=}
[ -n "$A" ] && [ -n "$B" ] && [ -d "$LIB_A" ] && [ -d "$LIB_B" ] || {
  echo "usage: sh bench_blp.sh <tagA>=<libA> <tagB>=<libB> [<configs>]"; exit 1; }
OUT=${OUT:-$D/output/blp-$(date +%Y%m%d-%H%M%S)}; mkdir -p "$OUT"
LOG=$OUT/run.log
ROUNDS=${ROUNDS:-3}
QUIET_CPU=${QUIET_CPU:-300}
CONFIGS=${3:-}
if [ -z "$CONFIGS" ]; then
  CONFIGS=$OUT/configs.txt
  cat > "$CONFIGS" <<'EOF'
model=mnl regime=claims U=15000 methods=simple,anderson,squarem threads=1,11 reps=2
model=mnl regime=census U=200000 methods=simple,anderson,squarem threads=1,11 reps=2
model=mnl regime=small U=2500000 methods=simple,anderson,squarem threads=1,11 reps=2
model=nl regime=claims U=15000 methods=simple,anderson,squarem threads=1,11 reps=2
model=nl regime=census U=200000 methods=simple,anderson,squarem threads=1,11 reps=2
model=nl regime=small U=2500000 methods=simple,anderson,squarem threads=1,11 reps=2
model=mxl regime=claims U=1500 methods=simple,anderson,squarem threads=1,11 reps=2
model=mxl regime=census U=20000 methods=simple,anderson,squarem threads=1,11 reps=2
model=mxl regime=small U=250000 methods=simple,anderson,squarem threads=1,11 reps=2
EOF
fi
ps_top() { ps -A -o %cpu=,comm= | sort -rn | head -4 | awk '{n=split($2,p,"/"); printf "%s%% %s; ", $1, p[n]}'; }
cpu_all() { ps -A -o %cpu= | awk '{s+=$1} END {printf "%d", s}'; }
load_avg() { if [ -r /proc/loadavg ]; then cut -d' ' -f1-3 /proc/loadavg; else sysctl -n vm.loadavg; fi; }
wait_quiet() {
  n=0
  while [ "$(cpu_all)" -gt "$QUIET_CPU" ] && [ $n -lt 60 ]; do sleep 10; n=$((n+1)); done
  [ $n -ge 60 ] && echo "  BUSY: still $(cpu_all)% after 10 min ($(ps_top))" >> "$LOG"
  return 0
}
run_one() {  # tag lib line_no round args
  wait_quiet
  echo "$(date +%H:%M:%S) $1 line $3 round $4 | load $(load_avg) | $(ps_top)" >> "$LOG"
  R_LIBS=$2 OMP_NUM_THREADS=11 Rscript "$D/bench_blp.R" tag=$1 $5 \
    out="$OUT/$1_$3_$4.rds" >> "$LOG" 2>&1 || echo "  FAILED: $1 line $3 round $4" >> "$LOG"
}
i=0
grep -v '^#' "$CONFIGS" | grep -v '^ *$' | while IFS= read -r line; do
  i=$((i + 1))
  r=1
  while [ $r -le "$ROUNDS" ]; do
    if [ $((r % 2)) -eq 1 ]; then
      run_one "$A" "$LIB_A" $i $r "$line"; run_one "$B" "$LIB_B" $i $r "$line"
    else
      run_one "$B" "$LIB_B" $i $r "$line"; run_one "$A" "$LIB_A" $i $r "$line"
    fi
    r=$((r + 1))
  done
done
echo "done: $OUT" | tee -a "$LOG"
