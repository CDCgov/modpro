#!/usr/bin/env python3

import argparse
import csv
import shutil
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--original", required=True, help="Raw HA1 PDB")
    p.add_argument("--relaxed", required=True, help="Rosetta relaxed PDB")
    p.add_argument("--original-score", required=True, help="Raw HA1 Rosetta scorefile")
    p.add_argument("--relaxed-score", required=True, help="Relaxed Rosetta scorefile")
    p.add_argument("--selected", required=True, help="Output selected PDB")
    p.add_argument("--out", required=True, help="Decision TSV")
    p.add_argument("--log", required=False, help="Optional results log file")
    return p.parse_args()


def read_lowest_rosetta_score(scorefile):
    scorefile = Path(scorefile)
    if not scorefile.exists():
        raise FileNotFoundError(f"Missing scorefile: {scorefile}")

    header = None
    best = None

    with scorefile.open() as handle:
        for line in handle:
            line = line.strip()
            if not line.startswith("SCORE:"):
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            if parts[1] in {"score", "total_score"}:
                header = parts[1:]
                continue

            if header is None:
                continue

            row = dict(zip(header, parts[1:]))

            value = None
            for key in ["total_score", "score"]:
                if key in row:
                    try:
                        value = float(row[key])
                        break
                    except ValueError:
                        pass

            if value is None:
                continue

            if best is None or value < best:
                best = value

    if best is None:
        raise ValueError(f"No usable Rosetta score found in {scorefile}")

    return best


def main():
    args = parse_args()

    original = Path(args.original)
    relaxed = Path(args.relaxed)
    selected = Path(args.selected)
    out = Path(args.out)

    if not original.exists():
        raise FileNotFoundError(f"Missing original PDB: {original}")
    if not relaxed.exists():
        raise FileNotFoundError(f"Missing relaxed PDB: {relaxed}")

    original_score = read_lowest_rosetta_score(args.original_score)
    relaxed_score = read_lowest_rosetta_score(args.relaxed_score)

    if relaxed_score < original_score:
        chosen = "relaxed"
        chosen_pdb = relaxed
        chosen_score = relaxed_score
    else:
        chosen = "original"
        chosen_pdb = original
        chosen_score = original_score

    selected.parent.mkdir(parents=True, exist_ok=True)
    out.parent.mkdir(parents=True, exist_ok=True)

    shutil.copyfile(chosen_pdb, selected)

    with out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "chosen",
            "chosen_pdb",
            "selected_pdb",
            "chosen_score",
            "original_pdb",
            "original_score",
            "relaxed_pdb",
            "relaxed_score",
        ])
        writer.writerow([
            chosen,
            chosen_pdb,
            selected,
            chosen_score,
            original,
            original_score,
            relaxed,
            relaxed_score,
        ])

    summary = "\n".join([
    f"chosen={chosen}",
    f"chosen_score={chosen_score}",
    f"selected_pdb={selected}",
    f"original_pdb={original}",
    f"original_score={original_score}",
    f"relaxed_pdb={relaxed}",
    f"relaxed_score={relaxed_score}",
    ]) + "\n"

    print(summary, end="")

    if args.log:
        log = Path(args.log)
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(summary)

if __name__ == "__main__":
    main()

