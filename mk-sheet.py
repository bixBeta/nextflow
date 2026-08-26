#!/usr/bin/env python3
"""
mk-sheet.py  —  Generate a sample-sheet.csv for the TREx RNA-seq Nextflow pipeline.

Usage:
    python mk-sheet.py <path>              # absolute paths in sheet (default)
    python mk-sheet.py <path> --names      # filenames only (use with --fastqs param)
    python mk-sheet.py <path> --se         # single-end (label, fastq1 only)
    python mk-sheet.py <path> -o out.csv   # custom output filename

<path> can be:
    - A Project_* directory  (e.g. Unaligned/Project_10488629)
    - An Unaligned directory (script finds Project_* inside)
    - Any directory containing Sample_* subdirectories
"""

import argparse
import csv
import re
import sys
from pathlib import Path


def derive_label(sample_dir_name: str) -> str:
    """
    Strip 'Sample_' prefix, then take all underscore-delimited tokens up to
    (but not including) the first purely numeric token.

    Example:
        Sample_8027D_BH10_10_10488629_23L335LT3_L7  →  8027D_BH10
        Sample_CTL_1_10488629_23L335LT3_L7           →  CTL
    """
    name = re.sub(r"^Sample_", "", sample_dir_name)
    tokens = name.split("_")
    label_tokens = []
    for tok in tokens:
        if tok.isdigit():
            break
        label_tokens.append(tok)
    return "_".join(label_tokens) if label_tokens else name


def find_fastqs(sample_dir: Path):
    """Return (R1_files, R2_files) lists sorted, searching recursively."""
    all_fqs = sorted(sample_dir.rglob("*.fastq.gz")) + sorted(sample_dir.rglob("*.fq.gz"))
    r1 = [f for f in all_fqs if re.search(r"_R1[_.]|_1\.f(ast)?q", f.name)]
    r2 = [f for f in all_fqs if re.search(r"_R2[_.]|_2\.f(ast)?q", f.name)]
    return r1, r2


def find_project_dir(root: Path) -> Path:
    """Resolve the directory that directly contains Sample_* subdirs."""
    if any(c.is_dir() and c.name.startswith("Sample_") for c in root.iterdir()):
        return root
    projects = [c for c in root.iterdir() if c.is_dir() and c.name.startswith("Project_")]
    if len(projects) == 1:
        return projects[0]
    if len(projects) > 1:
        print(f"Multiple Project_* dirs found — please point directly to one:", file=sys.stderr)
        for p in projects:
            print(f"  {p}", file=sys.stderr)
        sys.exit(1)
    print(f"No Sample_* or Project_* directories found under {root}", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Generate sample-sheet.csv for the TREx RNA-seq Nextflow pipeline."
    )
    parser.add_argument("path", help="Path to Project_* or Unaligned directory")
    parser.add_argument(
        "--names", action="store_true",
        help="Write filenames only (use with --fastqs in Nextflow). Default: absolute paths."
    )
    parser.add_argument(
        "--se", action="store_true",
        help="Single-end mode: output label + fastq1 only (no fastq2 column)."
    )
    parser.add_argument(
        "-o", "--output", default="sample-sheet.csv",
        help="Output CSV filename (default: sample-sheet.csv)"
    )
    args = parser.parse_args()

    root = Path(args.path).resolve()
    if not root.exists():
        print(f"Path not found: {root}", file=sys.stderr)
        sys.exit(1)

    project_dir = find_project_dir(root)
    sample_dirs = sorted(
        d for d in project_dir.iterdir()
        if d.is_dir() and d.name.startswith("Sample_")
    )

    if not sample_dirs:
        print(f"No Sample_* directories found in {project_dir}", file=sys.stderr)
        sys.exit(1)

    rows = []
    warnings = []

    for sd in sample_dirs:
        label = derive_label(sd.name)
        r1_files, r2_files = find_fastqs(sd)

        if not r1_files:
            warnings.append(f"  WARNING: no R1 fastq found in {sd.name} — skipped")
            continue

        if len(r1_files) > 1:
            warnings.append(
                f"  WARNING: multiple R1 files in {sd.name} — using first: {r1_files[0].name}"
            )

        r1 = r1_files[0]
        fastq1 = r1.name if args.names else str(r1)

        if args.se:
            rows.append({"label": label, "fastq1": fastq1})
        else:
            if not r2_files:
                warnings.append(f"  WARNING: no R2 fastq found in {sd.name} — skipped")
                continue
            r2 = r2_files[0]
            fastq2 = r2.name if args.names else str(r2)
            rows.append({"label": label, "fastq1": fastq1, "fastq2": fastq2})

    if warnings:
        print("\n".join(warnings), file=sys.stderr)

    if not rows:
        print("No samples written — check warnings above.", file=sys.stderr)
        sys.exit(1)

    fieldnames = ["label", "fastq1"] if args.se else ["label", "fastq1", "fastq2"]
    out_path = Path(args.output)
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    mode_label = "SE" if args.se else "PE"
    path_label = "filenames only (--fastqs mode)" if args.names else "absolute paths"
    print(f"\n{len(rows)} samples written to {out_path}  [{mode_label}, {path_label}]\n")

    # Preview first 3 rows
    print(f"  {'label':<25} {'fastq1':<50}" + ("" if args.se else " fastq2"))
    print(f"  {'-'*25} {'-'*50}" + ("" if args.se else f" {'-'*30}"))
    for row in rows[:3]:
        fq1 = Path(row["fastq1"]).name if not args.names else row["fastq1"]
        line = f"  {row['label']:<25} {fq1:<50}"
        if not args.se:
            fq2 = Path(row["fastq2"]).name if not args.names else row["fastq2"]
            line += f" {fq2}"
        print(line)
    if len(rows) > 3:
        print(f"  ... ({len(rows) - 3} more)")
    print()


if __name__ == "__main__":
    main()
