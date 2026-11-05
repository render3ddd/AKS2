#!/usr/bin/env bash
# Time the corrected C++ (render3d/AKS, tag v1.1) on the RUST_PLAN.md §2 inputs.
#
#   tools/cpp_baselines.sh ../AKS            # path to a built AKS checkout (`make`)
#
# Prints a Markdown table to paste into RUST_PLAN.md §2. Runs each case
# $REPS times (default 3) and reports the fastest wall-clock time. The C++
# uses every hardware thread for the Zn[x] congruence loop; close other apps.
set -euo pipefail

AKS=${1:?usage: $0 path/to/AKS}
REPS=${REPS:-3}
ZNX="$AKS/devVR/LenstraZnx.out"
Z="$AKS/devVR/LenstraZ.out"
[[ -x $ZNX && -x $Z ]] || { echo "build the C++ first: (cd $AKS && make)" >&2; exit 2; }

best() {                      # best <binary> <n>  -> fastest seconds of $REPS runs
    local t min=
    for _ in $(seq "$REPS"); do
        t=$( { TIMEFORMAT=%R; time (cd "$AKS" && "$1" "$2" >/dev/null); } 2>&1 )
        if [[ -z $min ]] || awk "BEGIN{exit !($t < $min)}"; then min=$t; fi
    done
    echo "$min"
}

echo "Machine: $(uname -sm), $(getconf _NPROCESSORS_ONLN) hardware threads"
echo
echo "| n | Zn[x] (NTL) s | Z (binary segmentation) s |"
echo "|---|---|---|"
for n in 23456789 1000000007; do
    echo "| $n | $(best "$ZNX" "$n") | $(best "$Z" "$n") |"
done
# The 16-digit case takes ~100 s per Zn[x] run; Z is skipped (far slower).
if [[ ${LONG:-0} == 1 ]]; then
    echo "| 1003026954441971 | $(REPS=1 best "$ZNX" 1003026954441971) | — |"
fi
