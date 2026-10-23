#!/bin/sh
# A sweep of the rows a thread takes per dispatch in the multinomial and nested
# logits' loops over choice situations (CHOICER_CHUNK_ROWS, src/choicer_derivs.h;
# 256 by default).
#
#   sh _benchmarks/scale/chunk_sweep.sh <choicer source tree>
#   env: OUT (default _benchmarks/scale/output/chunk-<date>), ROWS (default
#   "1 64 256 1024"), THREADS (default "1,5,11"), ROUNDS (default 4), CASES
#   (default all five, as model:regime:units:nests, units as bench_mnl_nl.R's U)
#
# Those loops hand a thread about CHOICER_CHUNK_ROWS stacked rows per dispatch,
# at least one situation, and fewer where the largest situation has more than
# eight times the mean rows (ROWS=1 builds the schedule of one situation per
# dispatch that 0.3.0 replaced). The value is fixed when the package is
# compiled, so the script builds one library per value from a copy of the tree,
# with the machine's Makevars (tools::makevars_user()) plus the define through
# R_MAKEVARS_USER, checks that the tree uses the define, finds it in each
# build's log, and checks that each library loads from where it was installed
# and has OpenMP (and the summary, that every run loaded its own). It times
# bench_mnl_nl.R's gradient, BHHH, Hessian and predictions on
# the same synthetic data with each: per round, every case under every library,
# the libraries in reverse order in alternate rounds, one R process per run.
# Cases: the multinomial logit on 2.5 million four-alternative situations
# (small) and on the census- and claims-style designs; the nested logit on the
# small design (two nests) and the census-style one (50 nests). Results: one
# .rds per run, run.log (with the load before each run), and summary.txt: per
# case, kernel and thread count, each library's median time and its ratio to
# 256's in every round. A run resumes in an OUT that holds results, skipping the
# runs it has. One thread's results do not depend on the value; the package's
# tests assume the default. Decision rule: README.md.
D=$(cd "$(dirname "$0")" && pwd)
SRC=$1
[ -n "$SRC" ] && [ -f "$SRC/DESCRIPTION" ] || { echo "usage: sh chunk_sweep.sh <choicer source tree>"; exit 1; }
grep -q CHOICER_CHUNK_ROWS "$SRC/src/choicer_derivs.h" ||
  { echo "$SRC does not use CHOICER_CHUNK_ROWS (a tree from before 0.3.0?)"; exit 1; }
OUT=${OUT:-$D/output/chunk-$(date +%Y%m%d-%H%M%S)}; mkdir -p "$OUT"
ROWS=${ROWS:-"1 64 256 1024"}; THREADS=${THREADS:-1,5,11}; ROUNDS=${ROUNDS:-4}
case " $ROWS " in *" 256 "*) ;; *) echo "ROWS must include 256, the default the summary compares with"; exit 1 ;; esac
LOG=$OUT/run.log
CASES=${CASES:-"mnl:small:2500000:2 mnl:census:20000:2 mnl:claims:15000:2 nl:small:2500000:2 nl:census:20000:50"}
load_avg() { if [ -r /proc/loadavg ]; then cut -d' ' -f1-3 /proc/loadavg; else sysctl -n vm.loadavg; fi; }

# One library per value, each from a fresh copy of the tree (make keeps object
# files whose sources did not change, whatever the flags)
MK_USER=$(Rscript --no-echo -e 'cat(tools::makevars_user())')
for v in $ROWS; do
  lib=$OUT/lib-$v
  if [ ! -d "$lib/choicer" ]; then
    s=$OUT/src-$v; rm -rf "$s" "$lib"; mkdir -p "$s" "$lib"
    for x in DESCRIPTION NAMESPACE R src inst data man; do [ -e "$SRC/$x" ] && cp -R "$SRC/$x" "$s/"; done
    rm -f "$s"/src/*.o "$s"/src/*.so "$s"/src/*.dll
    { [ -n "$MK_USER" ] && cat "$MK_USER"; echo "CPPFLAGS += -DCHOICER_CHUNK_ROWS=$v"; } > "$OUT/Makevars-$v"
    R_MAKEVARS_USER="$OUT/Makevars-$v" R CMD INSTALL --no-docs --no-multiarch --preclean \
      --library="$lib" "$s" > "$OUT/install-$v.log" 2>&1 ||
      { echo "building lib-$v failed; see $OUT/install-$v.log"; exit 1; }
    rm -rf "$s"
  fi
  grep -q -- "-DCHOICER_CHUNK_ROWS=$v" "$OUT/install-$v.log" ||
    { echo "lib-$v: the define is not in $OUT/install-$v.log"; exit 1; }
  # the library to load comes as an argument, not from R_LIBS, which a user
  # Renviron can override (the benchmarks would then all time that library)
  R_LIBS="$lib" Rscript --no-echo -e 'lib <- commandArgs(TRUE)[1]
    suppressMessages(library(choicer, lib.loc = lib))
    stopifnot(normalizePath(strsplit(Sys.getenv("R_LIBS"), ":")[[1]][1]) == normalizePath(lib),
              normalizePath(dirname(getNamespaceInfo("choicer", "path"))) == normalizePath(lib),
              isTRUE(thread_info()$openmp_enabled))' "$lib" ||
    { echo "lib-$v does not load from $lib with OpenMP (does a user Renviron set R_LIBS?)"; exit 1; }
done

echo "==== start $(date) | rows $ROWS | threads $THREADS | rounds $ROUNDS | source $SRC | load $(load_avg)" >> "$LOG"
r=1
while [ "$r" -le "$ROUNDS" ]; do
  order=$ROWS
  [ $((r % 2)) = 0 ] && order=$(echo "$ROWS" | awk '{for (i = NF; i >= 1; i--) printf "%s ", $i}')
  for c in $CASES; do
    set -- $(echo "$c" | tr ':' ' ')
    for v in $order; do
      f=$OUT/$1_$2_ch${v}_r$r.rds
      [ -f "$f" ] && continue
      echo "---- $(date +%H:%M:%S) ch$v $1 $2 r$r | load $(load_avg)" >> "$LOG"
      env R_LIBS="$OUT/lib-$v" OPENBLAS_NUM_THREADS=1 CHOICER_SCALE_OUT="$OUT" \
        Rscript --no-echo --no-restore "$D/bench_mnl_nl.R" tag=ch$v model=$1 regime=$2 U=$3 \
        nests=$4 threads=$THREADS reps=3 what=grad,bhhh,hess,predict out="$f" >> "$LOG" 2>&1
    done
  done
  r=$((r + 1))
done
echo "==== done $(date)" >> "$LOG"

# The summary, from a file: R's front end passes -e code through echo, which
# eats backslashes on macOS's sh and dash
cat > "$OUT/summary.R" <<'EOR'
out <- commandArgs(TRUE)[1]
fs <- list.files(out, "_r[0-9]+[.]rds$", full.names = TRUE)
rows <- do.call(rbind, lapply(fs, function(f) {
  x <- readRDS(f); i <- attr(x, "info")
  # each run must have loaded its own library
  if (basename(i$lib) != paste0("lib-", sub("ch", "", i$tag)))
    stop(basename(f), " loaded ", i$lib)
  x <- x[!is.na(x$sec), ]
  if (!nrow(x)) return(NULL)  # every kernel failed
  a <- aggregate(sec ~ what + threads, data = x, FUN = median)
  cbind(tag = i$tag, case = paste0(i$model, ":", i$regime),
        round = gsub("^.*_r|[.]rds$", "", basename(f)), a)
}))
tags <- unique(rows$tag[order(as.numeric(sub("ch", "", rows$tag)))])
for (g in split(rows, interaction(rows$case, rows$what, rows$threads, drop = TRUE))) {
  med <- sapply(tags, function(t) median(g$sec[g$tag == t]))
  rounds <- sort(unique(as.integer(g$round)))
  # ratios to 256 round by round; NA where either run is missing
  sec_of <- function(t, r) { v <- g$sec[g$tag == t & as.integer(g$round) == r]; if (length(v)) v else NA }
  rr <- sapply(setdiff(tags, "ch256"), function(t) {
    q <- vapply(rounds, function(r) sec_of(t, r) / sec_of("ch256", r), 0)
    sprintf("%s %s", t, paste(sprintf("%.2f", q), collapse = " "))
  })
  cat(sprintf("%-11s %-8s t=%-2d %s | vs ch256: %s\n", g$case[1], g$what[1], g$threads[1],
              paste(sprintf("%s %.4f", tags, med), collapse = "  "), paste(rr, collapse = " ; ")))
}
EOR
Rscript --no-echo "$OUT/summary.R" "$OUT" > "$OUT/summary.txt"
cat "$OUT/summary.txt"
