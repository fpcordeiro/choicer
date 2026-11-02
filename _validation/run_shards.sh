#!/bin/sh
# ============================================================================
# Runs one mc.R run as K shards in the background. First a warm pass (mc.R
# --warm, once: it runs the library's compile hooks, because parallel
# Rcpp::sourceCpp() calls into one cache race; it also checks the selection
# and names the run directory), then K shards started with nohup, each from
# its own generated script <run>/logs/shard_<k>_of_<K>.sh (with .log and
# .pid beside it), so that each can be stopped by its pid or by its
# script's name. A shard's script stops its R process when it is sent TERM.
# Prints the pids and the commands to watch and stop the shards.
#
# Usage (R_LIBS must name the build under test first; POSIX sh):
#   R_LIBS=<lib> _validation/run_shards.sh --family=F3 --shards=K [--threads=T]
#       [--run=TAG] [--out=DIR] [other mc.R arguments]
# Each shard runs at T OpenMP threads (default 1). --run defaults to a time
# stamp. Rscript is $RSCRIPT, else the one on the PATH.
# ============================================================================
set -eu

die() { printf 'run_shards.sh: %s\n' "$*" >&2; exit 1; }
# Its argument single-quoted for sh.
q() { printf "'%s'" "$(printf '%s' "$1" | sed "s/'/'\\\\''/g")"; }
is_count() { case $1 in ''|*[!0-9]*) return 1 ;; esac; [ "$1" -ge 1 ]; }

[ -n "${R_LIBS:-}" ] || die "R_LIBS is not set: name the build under test (R_LIBS=<lib> $0 ...)"
VDIR=$(cd "$(dirname "$0")" && pwd)
RSCRIPT=${RSCRIPT:-Rscript}

# Take --shards, --threads and --run; keep the other arguments for mc.R.
shards= threads=1 run= family=
n=$#
i=0
while [ "$i" -lt "$n" ]; do
  a=$1
  shift
  i=$((i + 1))
  case $a in
    --shards=*) shards=${a#--shards=}; continue ;;
    --threads=*) threads=${a#--threads=}; continue ;;
    --run=*) run=${a#--run=}; continue ;;
    --shard=*) die "--shard is set per shard: use --shards=K" ;;
    --warm|--check-seeds|--warm=*|--check-seeds=*) die "$a is not an argument of a run" ;;
    --family=?*) family=${a#--family=} ;;
  esac
  set -- "$@" "$a"
done
[ -n "$family" ] || die "--family=F is required (an id, a list a,b, or all)"
is_count "$shards" || die "--shards=K is required, K a positive whole number"
is_count "$threads" || die "--threads=T must be a positive whole number"
[ -n "$run" ] || run=$(date +%Y%m%d-%H%M%S)

# The warm pass: compile hooks, the selection checked, the run directory.
warm=$(OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 "$RSCRIPT" "$VDIR/mc.R" --warm --run="$run" "$@" 2>&1) || {
  printf '%s\n' "$warm" >&2
  die "the warm pass failed"
}
printf '%s\n' "$warm" | grep -v '^RUN_DIR=' || true
run_dir=$(printf '%s\n' "$warm" | sed -n 's/^RUN_DIR=//p' | tail -n 1)
[ -n "$run_dir" ] && [ -d "$run_dir" ] || die "the warm pass named no run directory"
logs=$run_dir/logs
mkdir -p "$logs"

args=
for a in "$@"; do args="$args $(q "$a")"; done

# No shard starts while any shard of this run still runs.
k=1
while [ "$k" -le "$shards" ]; do
  base=$logs/shard_${k}_of_${shards}
  [ ! -f "$base.pid" ] || ! kill -0 "$(cat "$base.pid")" 2>/dev/null ||
    die "shard $k/$shards of run $run is running already (pid $(cat "$base.pid")); no shard started"
  k=$((k + 1))
done

k=1
while [ "$k" -le "$shards" ]; do
  base=$logs/shard_${k}_of_${shards}
  {
    printf '#!/bin/sh\n'
    printf '# Shard %d/%d of run %s, written by run_shards.sh on %s.\n' "$k" "$shards" "$run" "$(date)"
    printf '# TERM stops its R process; the log is %s.log.\n' "$base"
    printf 'export R_LIBS=%s\n' "$(q "$R_LIBS")"
    printf 'export OMP_NUM_THREADS=%s OMP_THREAD_LIMIT=%s\n' "$threads" "$threads"
    printf '%s %s --run=%s --threads=%s --shard=%d/%d%s &\n' "$(q "$RSCRIPT")" "$(q "$VDIR/mc.R")" \
      "$(q "$run")" "$threads" "$k" "$shards" "$args"
    printf 'child=$!\n'
    printf "trap 'kill -TERM \$child 2>/dev/null' TERM\n"
    printf 'status=0\n'
    printf 'wait $child || status=$?\n'
    printf 'while kill -0 $child 2>/dev/null; do wait $child || status=$?; done\n'
    printf 'echo "shard %d/%d: exit status $status, $(date)"\n' "$k" "$shards"
    printf 'exit $status\n'
  } > "$base.sh"
  nohup sh "$base.sh" > "$base.log" 2>&1 &
  echo $! > "$base.pid"
  printf 'shard %d/%d: pid %s; script %s.sh; log %s.log\n' "$k" "$shards" "$!" "$base" "$base"
  k=$((k + 1))
done

# Stops every shard whose pid still runs its script (a pid can be reused).
{
  printf '#!/bin/sh\n'
  printf '# Stops the shards of run %s, written by run_shards.sh.\n' "$run"
  printf 'for f in %s/shard_*_of_%d.pid; do\n' "$(q "$logs")" "$shards"
  printf '  pid=$(cat "$f"); s=$(basename "$f" .pid).sh\n'
  printf '  if ps -p "$pid" -o command= 2>/dev/null | grep -q "$s"; then\n'
  printf '    kill -TERM "$pid" && echo "stopped $s (pid $pid)"\n'
  printf '  fi\n'
  printf 'done\n'
} > "$logs/stop_shards.sh"

cat <<EOF
run directory: $run_dir
watch:     tail -f $(q "$logs")/shard_*_of_${shards}.log
progress:  find $(q "$run_dir") -name 'rep_*.rds' | wc -l
stop one:  kill <pid>    (or: pkill -f $(q "$logs/shard_<k>_of_${shards}.sh"))
stop all:  sh $(q "$logs/stop_shards.sh")
resume:    rerun this command with --run=$(q "$run") (finished units are skipped)
report:    R_LIBS=$(q "$R_LIBS") Rscript $(q "$VDIR/report.R") --run=$(q "$run_dir")
EOF
