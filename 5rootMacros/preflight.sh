#!/bin/sh
# Checks the machine and the run folder before the analysis chain starts, and stops with one line per problem instead of failing halfway through.
# Run from the run folder by run.sh and plot.sh. PREFLIGHT=0 skips it; PREFLIGHT_DONE=1 (set here) keeps plot.sh from repeating it when run.sh called it.
# ACLIC= (empty) skips the compile test, like run.sh and plot.sh.
[ "$PREFLIGHT" = "0" ] && exit 0
[ -n "$PREFLIGHT_DONE" ] && exit 0

bad=0
fail() { echo "preflight: $*" >&2; bad=1; }

command -v root >/dev/null 2>&1 || fail "root is not in PATH (module load root?)"
if command -v root-config >/dev/null 2>&1; then
    feat=$(root-config --features 2>/dev/null)
    for f in dataframe imt; do
        case " $feat " in *" $f "*) ;; *) fail "this ROOT ($(root-config --version 2>/dev/null)) lacks the $f feature" ;; esac
    done
fi

# ACLiC compile test: catches a missing include path (the vdt headers on the Mac) in a few seconds instead of in the first macro
if [ "${ACLIC-+O}" != "" ] && command -v root >/dev/null 2>&1; then
    tmp=$(mktemp -d 2>/dev/null || echo /tmp/preflight.$$) ; mkdir -p "$tmp"
    printf '#include <ROOT/RDataFrame.hxx>\nvoid pf(){ ROOT::RDataFrame df(1); (void)df; }\n' > "$tmp/pf.cpp"
    (cd "$tmp" && root -l -b -q "pf.cpp+" > out.log 2>&1) || fail "ACLiC cannot compile a macro here: $(grep -m1 -iE 'error' "$tmp/out.log" | cut -c1-140) (ACLIC= runs interpreted)"
    rm -rf "$tmp"
fi

mkdir -p cuts plots 2>/dev/null
for d in cuts plots; do
    t="$d/.preflight.$$"; ( : > "$t" ) 2>/dev/null && rm -f "$t" || fail "$d/ is not writable"
done
ls output*.root >/dev/null 2>&1 || fail "no output*.root in $(pwd)"
[ -f settings.csv ] || fail "settings.csv is missing in $(pwd)"

# free disk: the caches and cuts are about as large as the sim output; ask for twice that
need_kb=$(du -sk output*.root 2>/dev/null | awk '{s+=$1} END{print s*2}')
free_kb=$(df -Pk . 2>/dev/null | awk 'NR==2{print $4}')
if [ -n "$need_kb" ] && [ -n "$free_kb" ] && [ "$free_kb" -lt "$need_kb" ]; then fail "free disk $((free_kb / 1024)) MB, need at least $((need_kb / 1024)) MB (2x output*.root)"; fi

ncpu=${NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null)}
mem_kb=$(awk '/MemTotal/{print $2}' /proc/meminfo 2>/dev/null)
[ -z "$mem_kb" ] && mem_kb=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1024 ))
echo "preflight: ${ncpu:-?} cores, $((mem_kb / 1048576)) GB memory, $((free_kb / 1048576)) GB free disk, ROOT $(root-config --version 2>/dev/null)"
if [ -n "$ncpu" ] && [ "$mem_kb" -gt 0 ] && [ "$mem_kb" -lt $((ncpu / 8 * 4 * 1048576)) ]; then
    echo "preflight: warning: under 4 GB of memory per 8 threads, a 100M run may be killed" >&2
fi

[ "$bad" = "0" ] || { echo "preflight: failed (PREFLIGHT=0 skips these checks)" >&2; exit 1; }
