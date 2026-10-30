#!/usr/bin/env python3
"""
Turn a Rosetta residue_energy_breakdown file into the parquet that feeds
protein_modeling.relaxed_rosetta_per_resi_energy.

Writes the parquet DIRECTLY. There is no TSV round trip, because the target
is parquet and routing typed data through a text file is what caused every
bug here on 2026-10-27:

  * The script printed diagnostics to stdout while the rule redirected
    stdout into the data file, AND wrote the same path itself via to_csv.
    The surviving first line was the input path, so the converter took that
    path as the single column name. Impala reported
    "incompatible Parquet schema for column resi1. Column type: BIGINT,
     Parquet schema: optional byte_array .../<hash>.res_breakdown.txt".
  * to_csv used header=False, so even the correct file had no column names.
  * Writing a header, then re-inferring it, is a coin flip: a headerless TSV
    fed to a header-inferring reader silently eats the first residue as the
    header and loses a row per variant.

Writing parquet directly removes all three. dtypes are set here and go
straight to disk, so nothing can be re-guessed.

A TSV is still written when --tsv is given, for humans. Nothing reads it.

Target schema, 15 columns. Impala resolves parquet columns by POSITION by
default, which is why the old pipeline's files worked despite being named
"0".."14", but the names are matched to the table as well so the data is
also correct under name resolution.

     #  Impala column        type     from Rosetta
     0  resi1                bigint   resi1
     1  pdbid1               string   pdbid1
     2  restype1             string   restype1
     3  resi2                string   resi2        "--" on one-body rows
     4  pdbid2               string   pdbid2       "--" on one-body rows
     5  restype2             string   restype2     "onebody" on one-body rows
     6  hbond_short_bb       double   hbond_sr_bb
     7  hbond_long_bb        double   hbond_lr_bb
     8  hbond_bb_sidechain   double   hbond_bb_sc
     9  hbond_sidechain      double   hbond_sc
    10  disulfide            double   dslf_fa13
    11  omega                double   omega
    12  total_energy         double   total
    13  variant_hash         string   -v
    14  version              string   -t

resi2 is a STRING and keeps Rosetta's "--". It is not an integer and must
not be coerced to one: one-body rows carry intrinsic per-residue energies
and have no second residue.
"""
import argparse
import sys

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# Rosetta column -> Impala column, in table order. Also defines what is kept.
COLUMN_MAP = [
    ("resi1", "resi1"),
    ("pdbid1", "pdbid1"),
    ("restype1", "restype1"),
    ("resi2", "resi2"),
    ("pdbid2", "pdbid2"),
    ("restype2", "restype2"),
    ("hbond_sr_bb", "hbond_short_bb"),
    ("hbond_lr_bb", "hbond_long_bb"),
    ("hbond_bb_sc", "hbond_bb_sidechain"),
    ("hbond_sc", "hbond_sidechain"),
    ("dslf_fa13", "disulfide"),
    ("omega", "omega"),
    ("total", "total_energy"),
]

STRING_COLUMNS = ["pdbid1", "restype1", "resi2", "pdbid2", "restype2",
                  "variant_hash", "version"]
DOUBLE_COLUMNS = ["hbond_short_bb", "hbond_long_bb", "hbond_bb_sidechain",
                  "hbond_sidechain", "disulfide", "omega", "total_energy"]

# Exactly the table's schema, so the file cannot drift from it.
ARROW_SCHEMA = pa.schema(
    [pa.field("resi1", pa.int64())]
    + [pa.field(n, pa.string()) for n in ("pdbid1", "restype1", "resi2",
                                          "pdbid2", "restype2")]
    + [pa.field(n, pa.float64()) for n in DOUBLE_COLUMNS]
    + [pa.field("variant_hash", pa.string()), pa.field("version", pa.string())]
)


def log(msg):
    """Diagnostics go to stderr, always. stdout is never data here."""
    print(msg, file=sys.stderr)


def build(input_file, variant_hash, version):
    raw = pd.read_table(input_file, engine="python", delimiter=r"\s+", header=0)
    log(f"read {len(raw)} rows, {len(raw.columns)} columns from {input_file}")

    missing = [src for src, _ in COLUMN_MAP if src not in raw.columns]
    if missing:
        raise SystemExit(
            f"ERROR: {input_file} is missing expected Rosetta column(s): {missing}\n"
            f"Columns present: {list(raw.columns)}"
        )

    df = pd.DataFrame()
    for src, dest in COLUMN_MAP:
        df[dest] = raw[src]
    df["variant_hash"] = variant_hash
    df["version"] = version

    # resi1 is the one integer. "--" would be a Rosetta change, not normal.
    df["resi1"] = pd.to_numeric(df["resi1"], errors="coerce").astype("Int64")
    if df["resi1"].isna().any():
        log(f"WARNING: resi1 has {int(df['resi1'].isna().sum())} non-numeric value(s)")

    for col in STRING_COLUMNS:
        df[col] = df[col].astype("string")
    for col in DOUBLE_COLUMNS:
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("float64")

    onebody = (df["restype2"] == "onebody").sum()
    log(f"rows: {len(df)} ({onebody} one-body, {len(df) - onebody} pairwise)")
    return df


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input_file", required=True,
                        help="Rosetta residue_energy_breakdown file")
    parser.add_argument("-v", "--variant_hash", required=True,
                        help="Alphanumeric STRING linked to a specific amino acid sequence")
    parser.add_argument("-t", "--timestamp", required=True, dest="version",
                        help="Run date, written to the 'version' column as a string")
    parser.add_argument("-o", "--output", required=True,
                        help="Output parquet path")
    parser.add_argument("--tsv", help="Optional TSV copy, for humans only")
    args = parser.parse_args(argv)

    df = build(args.input_file, args.variant_hash, args.version)

    table = pa.Table.from_pandas(df, schema=ARROW_SCHEMA, preserve_index=False)
    pq.write_table(table, args.output, version="1.0")
    log(f"wrote {args.output}")

    if args.tsv:
        df.to_csv(args.tsv, sep="\t", header=True, index=False)
        log(f"wrote {args.tsv}")

    for name, type_ in zip(table.schema.names, table.schema.types):
        log(f"    {str(type_):10} {name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

