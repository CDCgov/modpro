import argparse
from pathlib import Path

import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", required=True)
parser.add_argument("-v", "--variant_hash", required=False)
parser.add_argument("-t", "--timestamp_version", required=False)
parser.add_argument("-o", "--output_file", required=True)
args = parser.parse_args()


def prepare_contacts(input_file):
    data = pd.read_csv(
        input_file,
        low_memory=False,
        sep="\t",
        header=None,
    )
    return data.drop(columns=[0])


run_prep_df = prepare_contacts(args.input_file)

if args.variant_hash:
    run_prep_df[11] = args.variant_hash

if args.timestamp_version:
    run_prep_df[12] = pd.to_datetime(args.timestamp_version)

output_path = Path(args.output_file)
output_path.parent.mkdir(parents=True, exist_ok=True)

run_prep_df.to_csv(
    output_path,
    sep="\t",
    header=False,
    index=False,
)
