#!/usr/bin/env bash
# Thread scaling of the congruence tests (RUST_PLAN.md M4).
#
#   tools/scaling.sh [n] [max-threads]       # defaults: 23456789, all hardware threads
#
# Builds the release binary, then times `aks -j T n` for T = 1, 2, 4, ... and
# max-threads, best of $REPS runs (default 3). Prints a Markdown table.
#
# On Apple Silicon, macOS puts busy threads on the performance cores first,
# so the speed-up should be near-linear up to the P-core count (5 or 6 on an
# M3 Pro) and flatter after that, where the efficiency cores join in.
set -euo pipefail

N=${1:-23456789}
MAX=${2:-$(getconf _NPROCESSORS_ONLN)}
REPS=${REPS:-3}
cd "$(dirname "$0")/.."
cargo build --release --quiet
BIN=target/release/aks

best() {
    local t min=
    for _ in $(seq "$REPS"); do
        t=$("$BIN" -j "$1" "$N" | sed -n 's/^Time taken: \([0-9]*\) milliseconds$/\1/p')
        if [[ -z $min || $t -lt $min ]]; then min=$t; fi
    done
    echo "$min"
}

ts=(); t=1
while (( t < MAX )); do ts+=("$t"); t=$(( t * 2 )); done
ts+=("$MAX")

echo "n = $N on $(uname -sm), $MAX hardware threads"
echo
echo "| threads | time (ms) | speed-up | efficiency |"
echo "|---:|---:|---:|---:|"
base=
for t in "${ts[@]}"; do
    ms=$(best "$t")
    base=${base:-$ms}
    awk -v t="$t" -v ms="$ms" -v b="$base" \
        'BEGIN { s = b / ms; printf "| %d | %d | %.2f× | %.0f%% |\n", t, ms, s, 100 * s / t }'
done
