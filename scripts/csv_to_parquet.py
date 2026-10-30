#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


parser = argparse.ArgumentParser(
    description="Convert a tab-delimited file to parquet."
)

parser.add_argument(
    "input",
    help="Input .tsv, .txt, or .csv file"
)

parser.add_argument(
    "-o",
    "--output",
    help="Exact output parquet path"
)

args = parser.parse_args()

input_path = Path(args.input)

if not input_path.exists():
    raise FileNotFoundError(f"Input file not found: {input_path}")

if args.output:
    output_path = Path(args.output)
else:
    output_path = input_path.with_suffix(".parquet")

output_path.parent.mkdir(parents=True, exist_ok=True)

data = pd.read_csv(
    input_path,
    sep="\t"
)

data_pa_table = pa.Table.from_pandas(
    data,
    preserve_index=False
)

pq.write_table(
    data_pa_table,
    output_path,
    version="1.0"
)

print(f"Wrote parquet: {output_path}")
