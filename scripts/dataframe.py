#!/usr/bin/env python3
import pandas as pd
import sys
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", required=True,
        help="Input PDB file")
parser.add_argument("-v", "--variant_hash", required=False,
        help="Alphanumeric character STRING linked to a specific amino acid sequence")
parser.add_argument("-r", "--reference", required=False,
        help="Alphanumeric character STRING linked to a specific amino acid sequence")
parser.add_argument("-t", "--timestamp", required=False,
        help="System TIMESTAMP applied to distinct analytical applications")
args = parser.parse_args()
print(args.input_file)
input=(args.input_file)

def prepare_(test_file):
    with open(test_file) as parsed_file:
        parsed_file = pd.read_table(parsed_file,engine="python",delimiter="\s+",header=0)
        print(parsed_file)
        rosetta_df = pd.DataFrame(data=parsed_file).drop(columns=["score","time","description","p_aa_pp","linear_chainbreak","lk_ball_wtd","overlap_chainbreak","pro_close","rama_prepro","fa_dun"])
        #rosetta_df.columns=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
        rosetta_df.info(verbose=True)
        print(input)
    return rosetta_df

run_prep_df = prepare_(input)
variant_hash = args.variant_hash
reference = args.reference
timestamp  = args.timestamp
reference = args.reference

if reference:
    run_prep_df[17] = reference
if variant_hash:
    run_prep_df[18] = variant_hash
if timestamp:
    run_prep_df[19] = pd.to_datetime(timestamp)
pdb_path = os.path.dirname(args.input_file)
out_file=f"{pdb_path}/{variant_hash}"+".features.tsv"
run_prep_df.to_csv(out_file,sep ='\t',header=False,index=False)
print(run_prep_df)