#!/bin/bash
# Extract real putative DNMs (and optionally synthetic DNMs) from an All of Us
# trio VCF exported from the VDS. No SLURM on the AoU workspace: regions run
# in parallel on the local VM with xargs.
#
# Usage:
#   bash run_aou_extract.sh <vcf.bgz> <ped> <chrom> <outdir> [nproc] [region_size] [swapped_ped]
#
# The VCF must be local (gsutil cp it and its .tbi first) and tabix-indexed.
# Outputs <outdir>/<chrom>_<start>-<end>_real.tsv (+ _synthetic.tsv when a
# swapped pedigree is given), the layout preprocess_features.py expects.

set -euo pipefail

VCF=${1:?vcf}
PED=${2:?ped}
CHROM=${3:?chrom}
OUTDIR=${4:?outdir}
NPROC=${5:-4}
SIZE=${6:-20000000}
SWAPPED_PED=${7:-}

HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUTDIR/logs"

REGIONS="$OUTDIR/regions_${CHROM}.txt"
python "$HERE/scripts/generate_regions.py" --vcf "$VCF" --chroms "$CHROM" \
    --region-size "$SIZE" --output "$REGIONS"
echo "$(wc -l < "$REGIONS") regions for $CHROM"

run_region() {
    local region=$1 mode=$2
    local safe=${region/:/_}
    local out="$OUTDIR/${safe}_${mode}.tsv"
    [[ -s "$out" ]] && { echo "skip $out"; return 0; }
    local extra=(--real_dnms)
    [[ "$mode" == synthetic ]] && extra=(--swapped_ped "$SWAPPED_PED")
    python "$HERE/scripts/extract_dnm_features.py" --vcf "$VCF" --ped "$PED" --aou \
        "${extra[@]}" --region "$region" --output "$out.partial" \
        2> "$OUTDIR/logs/${safe}_${mode}.log"
    mv "$out.partial" "$out"
    echo "done $out"
}
export -f run_region
export HERE VCF PED OUTDIR SWAPPED_PED

MODES=(real)
[[ -n "$SWAPPED_PED" ]] && MODES+=(synthetic)
for mode in "${MODES[@]}"; do
    xargs -P "$NPROC" -I{} bash -c "run_region {} $mode" < "$REGIONS"
done
