import argparse
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from Bio.PDB import PDBParser, DSSP

pd.options.mode.chained_assignment = None

# Bio.PDB.DSSP tuple positions. Named so they cannot be transposed again.
DSSP_SECONDARY_STRUCTURE = 2
DSSP_RELATIVE_ASA = 3

# Reference CA atoms, chain A. Column name -> residue number.
DISTANCE_TARGETS = {
    "distance_to_98": 98,
    "distance_to_190": 190,
    "distance_to_194": 194,
    "distance_to_153": 153,
}

COLUMNS = [
    "site", "chain_id", "node_id", "x_coord", "y_coord", "z_coord",
    "residue_type", "distance_to_98", "distance_to_190", "distance_to_194",
    "distance_to_153", "rsa", "asa", "variant_hash",
]

ARROW_SCHEMA = pa.schema([
    pa.field("site", pa.int64()),
    pa.field("chain_id", pa.string()),
    pa.field("node_id", pa.string()),
    pa.field("x_coord", pa.float64()),
    pa.field("y_coord", pa.float64()),
    pa.field("z_coord", pa.float64()),
    pa.field("residue_type", pa.string()),
    pa.field("distance_to_98", pa.float64()),
    pa.field("distance_to_190", pa.float64()),
    pa.field("distance_to_194", pa.float64()),
    pa.field("distance_to_153", pa.float64()),
    pa.field("rsa", pa.float64()),
    pa.field("asa", pa.string()),
    pa.field("variant_hash", pa.string()),
])


def log(msg):
    print(msg, file=sys.stderr)


def reference_atoms(model):
    """CA atoms to measure against. Missing ones give NaN distances."""
    found = {}
    for column, resnum in DISTANCE_TARGETS.items():
        try:
            found[column] = model["A"][resnum]["CA"]
        except KeyError:
            found[column] = None
            log(f"note: chain A residue {resnum} CA not found; {column} will be NULL")
    return found


def dssp_values(dssp, chain_id, res_id):
    """(relative ASA as float, secondary structure as str) for one residue."""
    key = (chain_id, res_id)
    if key not in dssp:
        return np.nan, None
    entry = dssp[key]
    ss = entry[DSSP_SECONDARY_STRUCTURE]
    rsa = pd.to_numeric(entry[DSSP_RELATIVE_ASA], errors="coerce")
    return rsa, (None if ss in (None, "") else str(ss))


def build(pdb_filename, variant_hash):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("X", pdb_filename)
    model = structure[0]
    dssp = DSSP(model, pdb_filename)

    refs = reference_atoms(model)
    rows = []

    for chain in model:
        for res in chain:
            if not res.id[0].isspace():      # skip hetero/water
                continue
            rsa, secondary_structure = dssp_values(dssp, chain.id, res.id)
            for atom in res:
                vector = atom.get_vector()
                row = {
                    "site": res.id[1],
                    "chain_id": chain.id,
                    "node_id": atom.get_name(),
                    "x_coord": atom.coord[0],
                    "y_coord": atom.coord[1],
                    "z_coord": atom.coord[2],
                    "residue_type": res.resname,
                    "rsa": rsa,
                    "asa": secondary_structure,
                    "variant_hash": variant_hash,
                }
                for column, ref_atom in refs.items():
                    row[column] = (
                        (vector - ref_atom.get_vector()).norm()
                        if ref_atom is not None else np.nan
                    )
                rows.append(row)

    if not rows:
        raise SystemExit(f"ERROR: no standard residues found in {pdb_filename}")

    df = pd.DataFrame(rows, columns=COLUMNS)
    df["site"] = pd.to_numeric(df["site"], errors="coerce").astype("Int64")
    for col in ("chain_id", "node_id", "residue_type", "asa", "variant_hash"):
        df[col] = df[col].astype("string")
    for col in ("x_coord", "y_coord", "z_coord", "rsa", *DISTANCE_TARGETS):
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("float64")

    log(f"{len(df)} atoms across {df['chain_id'].nunique()} chain(s), "
        f"{df['site'].nunique()} residues")
    log(f"rsa: {int(df['rsa'].notna().sum())} of {len(df)} non-null")
    log(f"asa (secondary structure) values: {sorted(set(df['asa'].dropna()))}")
    return df


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input_pdb_file", required=True, help="Input PDB file")
    parser.add_argument("-v", "--variant_hash", required=True,
                        help="Alphanumeric STRING linked to a specific amino acid sequence")
    parser.add_argument("-o", "--out_p", help="Output directory (legacy; parquet named from the hash)")
    parser.add_argument("--parquet", help="Exact output parquet path (overrides -o)")
    parser.add_argument("--tsv", help="Exact output TSV path (overrides -o)")
    args = parser.parse_args(argv)

    if not args.parquet and not args.out_p:
        parser.error("give --parquet, or -o for the legacy directory layout")

    parquet_path = args.parquet or f"{args.out_p}/{args.variant_hash}_epitope.parquet"
    tsv_path = args.tsv or (f"{args.out_p}/{args.variant_hash}_epitope.tsv" if args.out_p else None)

    df = build(args.input_pdb_file, args.variant_hash)

    table = pa.Table.from_pandas(df, schema=ARROW_SCHEMA, preserve_index=False)
    pq.write_table(table, parquet_path, version="1.0")
    log(f"wrote {parquet_path}")

    if tsv_path:
        df.to_csv(tsv_path, sep="\t", header=True, index=False)
        log(f"wrote {tsv_path}")

    for name, type_ in zip(table.schema.names, table.schema.types):
        log(f"    {str(type_):10} {name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

