#!/usr/bin/env python3

import argparse
import socket
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract per-residue pLDDT from PDB B-factor column."
    )
    p.add_argument("--pdb", required=True, help="Input PDB")
    p.add_argument("--out", required=True, help="Output parquet")
    p.add_argument("--variant", required=True, help="Variant hash")
    p.add_argument("--run", required=True, help="Run ID")
    return p.parse_args()


def residue_key(line):
    return {
        "chain": line[21].strip() or "_",
        "residue_name": line[17:20].strip(),
        "residue_number": int(line[22:26].strip()),
        "insertion_code": line[26].strip() or "",
    }


def main():
    args = parse_args()

    pdb = Path(args.pdb)
    out = Path(args.out)

    if not pdb.exists():
        raise FileNotFoundError(f"Missing PDB: {pdb}")

    rows = {}

    with pdb.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            try:
                key_data = residue_key(line)
                b_factor = float(line[60:66].strip())
            except ValueError:
                continue

            key = (
                key_data["chain"],
                key_data["residue_name"],
                key_data["residue_number"],
                key_data["insertion_code"],
            )

            if key not in rows:
                rows[key] = {
                    "run_id": args.run,
                    "variant_hash": args.variant,
                    **key_data,
                    "plddt_sum": 0.0,
                    "atom_count": 0,
                }

            rows[key]["plddt_sum"] += b_factor
            rows[key]["atom_count"] += 1

    records = []
    for row in rows.values():
        atom_count = row["atom_count"]
        mean_plddt = row["plddt_sum"] / atom_count if atom_count else None
        records.append(
            {
                "run_id": row["run_id"],
                "variant_hash": row["variant_hash"],
                "chain": row["chain"],
                "residue_name": row["residue_name"],
                "residue_number": row["residue_number"],
                "insertion_code": row["insertion_code"],
                "atom_count": atom_count,
                "plddt": mean_plddt,
            }
        )

    if not records:
        raise ValueError(f"No ATOM/HETATM records with B-factors found in {pdb}")

    df = pd.DataFrame(records).sort_values(
        ["chain", "residue_number", "insertion_code"]
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    
    # Write parquet in Impala-compatible format (v1.0)
    table = pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(
        table, 
        out, 
        version='1.0',
        coerce_timestamps='ms'
    )
    
    # Also write CSV for reference
    df.to_csv(out.with_suffix(".csv"), index=False)

    print(f"host={socket.gethostname()}")
    print(f"run={args.run}")
    print(f"variant={args.variant}")
    print(f"pdb={pdb}")
    print(f"out={out}")
    print(f"wrote_rows={len(df)}")
    print("plddt_extract_complete=true")


if __name__ == "__main__":
    main()
