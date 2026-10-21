#!/bin/sh
# macOS only: build the malloc-counting shim (memcount.c) and write the Rmc
# front end the memory runs use, both into output/.
#
#   sh _benchmarks/scale/make_rmc.sh
#   MC_SHIM=$PWD/_benchmarks/scale/output/memcount.dylib \
#     _benchmarks/scale/output/Rmc --no-echo --no-restore --file=<script> --args ...
#
# The shim is preloaded with DYLD_INSERT_LIBRARIES, which macOS strips from the
# environment of protected programs such as /bin/sh, which runs R's front end
# (bin/R). Rmc is a copy of this machine's front end that sets the variable
# from MC_SHIM just before it starts R's binary. Rewrite it after R is updated.
set -e
D=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$D/output"
[ "$(uname)" = Darwin ] || { echo "make_rmc.sh: the shim is for macOS"; exit 1; }
clang -O2 -dynamiclib -o "$D/output/memcount.dylib" "$D/memcount.c"
FRONT="$(R RHOME)/bin/R"
n=$(grep -c '^[[:space:]]*exec "\${R_binary}"' "$FRONT" || true)
[ "$n" -ge 1 ] || { echo "make_rmc.sh: no exec of \${R_binary} in $FRONT"; exit 1; }
awk '/^[[:space:]]*exec "\$\{R_binary\}"/ && !done {
       print "  if test -n \"${MC_SHIM}\"; then DYLD_INSERT_LIBRARIES=\"${MC_SHIM}\"; export DYLD_INSERT_LIBRARIES; fi"
       done = 1 }
     { print }' "$FRONT" > "$D/output/Rmc"
grep -q DYLD_INSERT_LIBRARIES "$D/output/Rmc" || { echo "make_rmc.sh: could not insert the shim line"; exit 1; }
chmod +x "$D/output/Rmc"
echo "wrote $D/output/memcount.dylib and $D/output/Rmc (from $FRONT)"
