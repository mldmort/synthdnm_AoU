#!/usr/bin/env python3
"""
Apply score cutoffs and region exclusions to classified candidate DNMs.

Reads one or more Parquets from classify.py, keeps real putative DNMs
(truth=0, or all rows if there is no truth column) that pass a per-type
score cutoff (SNV / indel, by indel_flag) and an optional cohort AC limit,
then removes calls whose anchor falls in the exclusion BEDs:

    genomicSuperDups, simpleRepeat, rmsk_simple_repeat   (resources/regions)

Same rule as the SPARK WGS 2026_08 release (James Guevara's
region_filter.py): each BED is merged (overlapping and touching intervals),
and a call is excluded if its 0-based anchor POS-1 lies in a merged
interval. Only the anchor is tested, not the full indel span. Counts per
BED are sequential, in the order above.

If the input also holds synthetic DNMs (truth=1), the fraction passing the
cutoff, and passing and lying outside all regions, is reported per type as
an upper bound on sensitivity. With --ped, calls per child are summarized
over every child in the pedigree (children with no calls count as zero).

Usage:
    python region_filter.py \\
        --input scored_odd_snv.parquet scored_odd_indel.parquet \\
        --snv_cutoff 0.3 --indel_cutoff 0.7 --max_ac 1 \\
        --ped trios.ped \\
        --output dnm_calls.parquet \\
        --excluded dnm_excluded_regions.parquet \\
        --summary dnm_calls_summary.json
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import polars as pl

REGION_NAMES = ["genomicSuperDups", "simpleRepeat", "rmsk_simple_repeat"]
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Apply score cutoffs and region exclusions to classified DNMs."
    )
    parser.add_argument("--input", required=True, nargs="+",
                        help="Parquet(s) from classify.py (e.g. one per SNV/indel model)")
    parser.add_argument("--output", required=True,
                        help="Output Parquet of passing calls outside all regions")
    parser.add_argument("--snv_cutoff", type=float, default=0.5,
                        help="Minimum synthdnm_prob for SNVs (default: 0.5)")
    parser.add_argument("--indel_cutoff", type=float, default=0.5,
                        help="Minimum synthdnm_prob for indels (default: 0.5)")
    parser.add_argument("--max_ac", type=int, default=None,
                        help="Keep only calls with cohort AC <= this (e.g. 1)")
    parser.add_argument("--regions", default=str(Path(__file__).resolve().parent.parent
                                                 / "resources" / "regions"),
                        help="Directory with <name>.bed[.gz] (default: repo resources/regions)")
    parser.add_argument("--ped", default=None,
                        help="PED file; summarize calls per child over all its children")
    parser.add_argument("--excluded", default=None,
                        help="Optional Parquet of passing calls removed by a region, with "
                             "one in_<region> flag per BED")
    parser.add_argument("--summary", default=None, help="Optional JSON summary")
    return parser.parse_args()


def find_bed(regions_dir, name):
    for suffix in (".bed.gz", ".bed"):
        path = Path(regions_dir) / f"{name}{suffix}"
        if path.exists():
            return path
    sys.exit(f"Missing {name}.bed[.gz] in {regions_dir}")


def load_merged(path):
    """Read a BED and return {chrom: (starts, ends)} with overlapping and
    touching intervals merged (start <= previous end), 0-based half-open."""
    bed = pl.read_csv(path, separator="\t", has_header=False, comment_prefix="#",
                      columns=[0, 1, 2], new_columns=["chrom", "start", "end"],
                      schema_overrides={"chrom": pl.Utf8, "start": pl.Int64, "end": pl.Int64})
    merged = {}
    for (chrom,), g in bed.sort(["chrom", "start"]).group_by(["chrom"], maintain_order=True):
        s = g["start"].to_numpy()
        e = g["end"].to_numpy()
        if np.any(s < 0) or np.any(e <= s):
            sys.exit(f"Invalid BED interval in {path} ({chrom})")
        new = np.r_[True, s[1:] > np.maximum.accumulate(e)[:-1]]
        merged[chrom] = (s[new], np.maximum.reduceat(e, np.flatnonzero(new)))
    return merged


def in_regions(chroms, positions, merged):
    """True where the 0-based anchor (POS-1) lies inside a merged interval."""
    pos = positions - 1
    flag = np.zeros(len(pos), dtype=bool)
    for chrom in np.unique(chroms):
        if chrom not in merged:
            continue
        idx = np.flatnonzero(chroms == chrom)
        starts, ends = merged[chrom]
        i = np.searchsorted(starts, pos[idx], side="right") - 1
        flag[idx] = (i >= 0) & (pos[idx] < ends[np.maximum(i, 0)])
    return flag


def describe_calls(df, children):
    """Counts, per-child distribution, indel fraction and SNV Ts/Tv."""
    out = {"calls": df.height}
    if df.height:
        snv = df.filter(pl.col("indel_flag") == 0)
        is_ts = np.array([(r, a) in TRANSITIONS
                          for r, a in zip(snv["REF"].to_list(), snv["ALT_specific"].to_list())],
                         dtype=bool)
        out["indel_fraction"] = round(float(df["indel_flag"].mean()), 4)
        out["snv_ts_tv"] = round(float(is_ts.sum() / max(1, (~is_ts).sum())), 3)
    if children is not None:
        per_child = (pl.DataFrame({"SAMPLE": children})
                     .join(df.group_by("SAMPLE").len(), on="SAMPLE", how="left")
                     .fill_null(0)["len"])
        out["children"] = len(children)
        out["per_child"] = {
            "mean": round(float(per_child.mean()), 3),
            "median": float(per_child.median()),
            "p99": float(per_child.quantile(0.99)),
            "max": int(per_child.max()),
            "pct_zero": round(100 * float((per_child == 0).mean()), 2),
        }
        out["snv_per_child"] = round((df["indel_flag"] == 0).sum() / len(children), 3)
        out["indel_per_child"] = round((df["indel_flag"] == 1).sum() / len(children), 3)
    return out


def main():
    args = parse_args()

    df = pl.concat([pl.read_parquet(p) for p in args.input], how="diagonal_relaxed")
    print(f"Read {df.height:,} rows from {len(args.input)} file(s)", file=sys.stderr)
    for col in ("#CHROM", "POS", "SAMPLE", "AC", "indel_flag", "synthdnm_prob"):
        if col not in df.columns:
            sys.exit(f"Missing column {col}; expected classify.py output")

    cutoff = (pl.when(pl.col("indel_flag") == 1).then(args.indel_cutoff)
              .otherwise(args.snv_cutoff))
    df = df.with_columns((pl.col("synthdnm_prob") >= cutoff).alias("passed"))

    regions = {}
    for name in REGION_NAMES:
        path = find_bed(args.regions, name)
        regions[name] = {"path": str(path),
                         "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                         "merged": load_merged(path)}
    chroms = df["#CHROM"].to_numpy()
    positions = df["POS"].to_numpy().astype(np.int64)
    df = df.with_columns([
        pl.Series(f"in_{name}", in_regions(chroms, positions, regions[name]["merged"]))
        for name in REGION_NAMES
    ])
    df = df.with_columns(pl.any_horizontal([f"in_{n}" for n in REGION_NAMES]).alias("in_any"))

    has_truth = "truth" in df.columns
    real = df.filter(pl.col("truth") == 0) if has_truth else df
    if args.max_ac is not None:
        real = real.filter(pl.col("AC") <= args.max_ac)
    n_candidates = real.height
    calls = real.filter(pl.col("passed"))

    print(f"Candidates{f' with AC <= {args.max_ac}' if args.max_ac is not None else ''}: "
          f"{n_candidates:,}", file=sys.stderr)
    print(f"Passing cutoffs (SNV >= {args.snv_cutoff}, indel >= {args.indel_cutoff}): "
          f"{calls.height:,}", file=sys.stderr)
    removed = {}
    remaining = calls
    for name in REGION_NAMES:
        hit = remaining.filter(pl.col(f"in_{name}"))
        removed[name] = hit.height
        remaining = remaining.filter(~pl.col(f"in_{name}"))
        print(f"  - {name}: {hit.height:,} removed -> {remaining.height:,} remaining",
              file=sys.stderr)

    children = None
    if args.ped:
        ped = pl.read_csv(args.ped, separator="\t", has_header=False, comment_prefix="#",
                          new_columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
                          schema_overrides={c: pl.Utf8 for c in
                                            ["FID", "IID", "PAT", "MAT", "SEX", "PHENO"]})
        children = ped.filter(pl.col("PAT") != "0")["IID"].unique().to_list()

    summary = {
        "inputs": args.input,
        "snv_cutoff": args.snv_cutoff,
        "indel_cutoff": args.indel_cutoff,
        "max_ac": args.max_ac,
        "candidates": n_candidates,
        "regions": {n: {"path": r["path"], "sha256": r["sha256"], "removed": removed[n]}
                    for n, r in regions.items()},
        "before_regions": describe_calls(calls, children),
        "after_regions": describe_calls(remaining, children),
    }

    if has_truth and df.filter(pl.col("truth") == 1).height:
        syn = df.filter(pl.col("truth") == 1)
        summary["synthetic"] = {}
        for flag, label in ((0, "snv"), (1, "indel")):
            g = syn.filter(pl.col("indel_flag") == flag)
            if g.height:
                summary["synthetic"][label] = {
                    "n": g.height,
                    "pass_cutoff": round(float(g["passed"].mean()), 4),
                    "pass_cutoff_outside_regions": round(float((g["passed"] & ~g["in_any"]).mean()), 4),
                    "in_any_region": round(float(g["in_any"].mean()), 4),
                }

    for stage in ("before_regions", "after_regions"):
        s = summary[stage]
        line = f"{stage}: {s['calls']:,} calls"
        if "per_child" in s:
            p = s["per_child"]
            line += (f", per child mean {p['mean']} (median {p['median']}, p99 {p['p99']}, "
                     f"max {p['max']}), SNV {s['snv_per_child']} + indel {s['indel_per_child']}")
        if "snv_ts_tv" in s:
            line += f", indel fraction {s['indel_fraction']}, SNV Ts/Tv {s['snv_ts_tv']}"
        print(line, file=sys.stderr)
    for label, s in summary.get("synthetic", {}).items():
        print(f"synthetic {label}: pass cutoff {s['pass_cutoff']}, pass and outside regions "
              f"{s['pass_cutoff_outside_regions']}, in any region {s['in_any_region']}",
              file=sys.stderr)

    flag_cols = [f"in_{n}" for n in REGION_NAMES] + ["in_any", "passed"]
    remaining.drop(flag_cols).write_parquet(args.output)
    print(f"Wrote {remaining.height:,} calls to {args.output}", file=sys.stderr)
    if args.excluded:
        calls.filter(pl.col("in_any")).drop(["in_any", "passed"]).write_parquet(args.excluded)
        print(f"Wrote {calls.filter(pl.col('in_any')).height:,} region-excluded calls to "
              f"{args.excluded}", file=sys.stderr)
    if args.summary:
        Path(args.summary).write_text(json.dumps(summary, indent=2) + "\n")
        print(f"Wrote summary to {args.summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
