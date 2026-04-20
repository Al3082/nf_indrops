#!/usr/bin/env bash
#
# One-time migration of existing output trees to the new stage-based naming.
#   v2    : Library_N       -> S{stage}_v2_{n}
#   v3NPB : NPB_1..4         -> S11_NPB_4..7
#
# Usage:
#   bash rename_legacy_libraries.sh <output_dir>            # dry run (default)
#   bash rename_legacy_libraries.sh <output_dir> --apply    # actually rename

set -euo pipefail

OUT="${1:-}"
MODE="${2:---dry-run}"

if [[ -z "$OUT" || ! -d "$OUT" ]]; then
    echo "Usage: $0 <output_dir> [--apply]" >&2
    exit 1
fi

declare -A MAP=(
    [Library_1]=S8_v2_1
    [Library_2]=S8_v2_2
    [Library_4]=S10_v2_1
    [Library_5]=S10_v2_2
    [Library_6]=S10_v2_3
    [Library_17]=S12_v2_1
    [Library_18]=S12_v2_2
    [Library_19]=S12_v2_3
    [Library_27]=S14_v2_1
    [Library_28]=S14_v2_2
    [Library_29]=S14_v2_3
    [Library_30]=S14_v2_4
    [Library_34]=S16_v2_1
    [Library_35]=S16_v2_2
    [Library_36]=S16_v2_3
    [Library_37]=S16_v2_4
    [Library_41]=S18_v2_1
    [Library_42]=S18_v2_2
    [Library_43]=S18_v2_3
    [Library_44]=S18_v2_4
    [Library_48]=S20_v2_1
    [Library_49]=S20_v2_2
    [Library_50]=S20_v2_3
    [Library_51]=S20_v2_4
    [Library_52]=S20_v2_5
    [Library_56]=S22_v2_1
    [Library_57]=S22_v2_2
    [Library_58]=S22_v2_3
    [Library_59]=S22_v2_4
    [Library_60]=S22_v2_5
    [Library_61]=S22_v2_6
    [NPB_1]=S11_NPB_4
    [NPB_2]=S11_NPB_5
    [NPB_3]=S11_NPB_6
    [NPB_4]=S11_NPB_7
)

APPLY=0
[[ "$MODE" == "--apply" ]] && APPLY=1

do_mv() {
    local src="$1" dst="$2"
    if (( APPLY )); then
        if [[ -e "$dst" ]]; then
            echo "SKIP (target exists): $src -> $dst" >&2
            return
        fi
        mv -v -- "$src" "$dst"
    else
        echo "DRY: $src -> $dst"
    fi
}

# Process longer names first (so Library_17 is handled before Library_1, etc.)
mapfile -t OLDS < <(
    for k in "${!MAP[@]}"; do printf '%d\t%s\n' "${#k}" "$k"; done \
    | sort -k1,1nr -k2,2 | cut -f2-
)

# --- Phase 1: rename top-level library directories ---------------------------
for old in "${OLDS[@]}"; do
    new="${MAP[$old]}"
    [[ -d "$OUT/$old" ]] && do_mv "$OUT/$old" "$OUT/$new"
done

# --- Phase 2: rename files whose basename is "<old>.*" or "<old>__*" ---------
# (covers processed_fastqs/, cutadapt_trim_stats/, qc/fastqc/, etc.)
for old in "${OLDS[@]}"; do
    new="${MAP[$old]}"
    while IFS= read -r -d '' f; do
        dir=$(dirname "$f")
        base=$(basename "$f")
        newbase="${new}${base#${old}}"
        do_mv "$f" "$dir/$newbase"
    done < <(
        find "$OUT" -depth -type f \
             \( -name "${old}.*" -o -name "${old}__*" \) -print0
    )
done

if (( ! APPLY )); then
    echo
    echo "Dry run complete. Re-run with --apply to perform the renames."
fi
