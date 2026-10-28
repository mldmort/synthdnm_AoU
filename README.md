# SynthDNM

Classifier for filtering de novo mutations (DNMs) from whole-genome sequencing trios.

Trains on synthetic DNMs (inherited variants via pedigree swapping) as proxies for real de novos, and putative DNMs (child het, parents hom ref) which are mostly artifacts. The classifier learns to distinguish high-quality real DNMs from noise.

## Pipeline

```
1. swap_pedigree.py        — Generate swapped pedigree
2. generate_regions.py     — Shard genome into regions for parallel extraction
3. extract_dnm_features.py — Extract genotype features (FORMAT + INFO fields)
4. preprocess_features.py  — Derive features, correct haploid_flag, output Parquet
5. train.py                — Train XGBoost classifier (3 feature sets)
6. classify.py             — Score candidate DNMs
```

## Feature Sets

Defined in `scripts/features.toml` with an extends chain:

| Feature Set | Features | Use Case |
|---|---|---|
| **universal** | 21 (FORMAT-only) | Any caller (DeepVariant, DRAGEN, etc.) |
| **gatk** | 28 (+INFO quality metrics) | GATK HaplotypeCaller |
| **ssc** | 29 (+VQSLOD) | SSC-specific GATK with VQSR |

## Label Semantics

- `synthdnm_prob` = P(truth=1) = probability variant is a real DNM
- **High prob** (close to 1) = looks like a real de novo = **keep**
- **Low prob** (close to 0) = looks like artifact = **discard**

Training labels:
- truth=1: Synthetic DNMs — inherited variants used as proxies for real DNMs (high-quality)
- truth=0: Putative DNMs — child het, parents hom ref (vast majority are artifacts)

## Validation

Evaluated against [denovo-db v1.6.1](https://denovo-db.gs.washington.edu/) SSC validated DNMs.
See `results/validation_summary.json` for full metrics.

| Model | Validated DNMs (n=3,006) | + AC<=1 filter (n=68,687) |
|---|---|---|
| universal | 96.8% recalled | 98.7% real |
| gatk | 95.5% recalled | 97.9% real |
| ssc | 96.5% recalled | 98.6% real |

Note: Additional downstream filters (gnomAD AF, segmental duplications, simple repeats) not yet applied.

## Usage

### Extract features (SLURM array jobs)
```bash
sbatch run_extract_real.sh
sbatch run_extract_synthetic.sh
```

### Preprocess + train
```bash
python scripts/preprocess_features.py \
    --feature_dir output/features/ \
    --output output/training_data.parquet \
    --gatk --sample_n 50000 \
    --psam resources/SSC.psam \
    --par resources/hg38_par.bed

python scripts/train.py \
    --input output/training_data.parquet \
    --output_dir output/models/gatk/ \
    --feature_set gatk \
    --config scripts/features.toml \
    --cv_folds 5
```

### Classify candidates
```bash
python scripts/classify.py \
    --input candidates.parquet \
    --model output/models/gatk/model.json \
    --output scored.tsv \
    --threshold 0.5
```

## All of Us

Input is a per-chromosome trio VCF exported from the AoU v9 VDS with Hail
(multiallelics split, FT-failed genotypes set to no-call). Expected fields:

| Level | Fields |
|---|---|
| FORMAT | GT, GQ, AD, DP (=sum AD), AB, RGQ, FT, PS, VR (1 = from `variant_data`, 0 = from a reference block) |
| INFO | AC, AN, AF (trio cohort), AS_VETS, VETS_MODEL, MULTIALLELIC |

Reference-block calls carry only GT and GQ (banded), so most parental AD/DP
are missing; `--aou` keeps them null rather than reading cyvcf2's missing
sentinels as depths. There is no PL, so use a no-PL feature set. Keep
AC/AN/AF out of the features: synthetic positives are AC=2 by construction.

```bash
# trio table with child/father/mother (+ sex) columns -> 6-column PED
python scripts/make_ped.py --trios trios.tsv --sex sex --output trios.ped
python scripts/swap_pedigree.py trios.ped -o trios.swapped.ped

# VCF must be local and tabix-indexed (gsutil cp the .bgz and .tbi first)
bash run_aou_extract.sh trios_chr22.vcf.bgz trios.ped chr22 output/features 4 20000000 trios.swapped.ped
```

`run_aou_extract.sh` runs `extract_dnm_features.py --aou` per 20 Mb region in
parallel (no SLURM) and writes `<region>_real.tsv` / `<region>_synthetic.tsv`;
omit the swapped PED to extract real putative DNMs only. A region whose output
already exists is skipped.

Then build the training table and train. `--max_real_ac 2` keeps only rare
putative DNMs as negatives, matching the AC=2 synthetic positives; without it
most negatives sit at sites carried by other families and the model learns
frequency. With `--sample_n`, both classes are capped at the smaller one.

Three AoU feature sets (see comments in `scripts/features.toml`):
`aou_core` (child GQ/DP/AB/RGQ, parental GQ, flags), `aou_core_vets` (+ VETS
site score) and `aou_core_vets_vr` (+ `n_parent_VR`, the number of parents
with a variant_data record at the site).

```bash
python scripts/preprocess_features.py --feature_dir output/features \
    --output output/training_data.parquet --sample_n 50000 --max_real_ac 2 --psam trios.ped
for fs in aou_core aou_core_vets aou_core_vets_vr; do
    python scripts/train.py --input output/training_data.parquet \
        --output_dir output/models/$fs --feature_set $fs --config scripts/features.toml
done
```

## Requirements

- Python 3.12+
- polars, pyarrow, xgboost, scikit-learn, numpy, cyvcf2 (XGBoost needs pyarrow to read polars tables)
