#!/usr/bin/env python3
"""
Build a 6-column PED file (FID IID PAT MAT SEX PHENO) from a trio table.

The extractor and swap_pedigree.py need one row per child plus one row per
parent, with parents' sex set (1=father, 2=mother). Cohorts like All of Us
provide trios but no family IDs, so FID is derived from the parent pair:
full siblings share a family, half-siblings do not.

Input: TSV/CSV with a header containing child, father and mother ID columns,
and optionally child sex (1/2, M/F, male/female; anything else -> 0).
IDs must match the VCF sample names exactly.

Usage:
    python make_ped.py --trios trios.tsv --output trios.ped
    python make_ped.py --trios trios.csv --child person_id --father father_id \\
        --mother mother_id --sex sex_at_birth --output trios.ped
"""

import argparse
import csv
import sys

SEX_CODES = {"1": "1", "m": "1", "male": "1", "2": "2", "f": "2", "female": "2"}


def parse_args():
    parser = argparse.ArgumentParser(description="Build a PED file from a trio table.")
    parser.add_argument("--trios", required=True, help="Trio table (TSV or CSV, with header)")
    parser.add_argument("--child", default="child", help="Child ID column (default: child)")
    parser.add_argument("--father", default="father", help="Father ID column (default: father)")
    parser.add_argument("--mother", default="mother", help="Mother ID column (default: mother)")
    parser.add_argument("--sex", default=None, help="Child sex column (optional)")
    parser.add_argument("--output", help="Output PED file (default: stdout)")
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.trios, newline="") as f:
        dialect = csv.Sniffer().sniff(f.readline(), delimiters="\t,")
        f.seek(0)
        rows = list(csv.DictReader(f, dialect=dialect))

    children = {}
    parents = {}
    for row in rows:
        child, father, mother = row[args.child].strip(), row[args.father].strip(), row[args.mother].strip()
        sex = SEX_CODES.get(row[args.sex].strip().lower(), "0") if args.sex else "0"
        if child in children and children[child][:2] != (father, mother):
            sys.exit(f"Child {child} listed with two different parent pairs")
        children[child] = (father, mother, sex)
        for pid, psex in ((father, "1"), (mother, "2")):
            if parents.setdefault(pid, psex) != psex:
                sys.exit(f"{pid} appears as both father and mother")

    overlap = set(children) & set(parents)
    if overlap:
        print(f"Warning: {len(overlap)} samples are both a child and a parent; "
              f"they are written as children only", file=sys.stderr)

    out = open(args.output, "w") if args.output else sys.stdout
    families = {}
    for child, (father, mother, sex) in children.items():
        fid = families.setdefault((father, mother), f"FAM{len(families) + 1:06d}")
        print(f"{fid}\t{child}\t{father}\t{mother}\t{sex}\t0", file=out)
    for (father, mother), fid in families.items():
        for pid, psex in ((father, "1"), (mother, "2")):
            if pid not in children:
                print(f"{fid}\t{pid}\t0\t0\t{psex}\t0", file=out)
    if args.output:
        out.close()

    print(f"{len(children)} children in {len(families)} families", file=sys.stderr)


if __name__ == "__main__":
    main()
