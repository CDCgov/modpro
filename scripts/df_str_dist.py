#!/usr/bin/env python3
import pandas as pd
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", required=True,
        help="Input PDB file")
args = parser.parse_args()
print(args.input_file)
input=(args.input_file)

def prepare_(test_file):
    with open(test_file) as parsed_file:
        parsed_file = pd.read_table(parsed_file,engine="python",delimiter="\s+",header=None)
        #print(parsed_file)
        df = pd.DataFrame(data=parsed_file)
        df2=df[[0,1,2,3,4,5,7,8,9]]
    return df2
    print(df2)

run_prep_df = prepare_(input)
#out_file=f"/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/cdp_archive/{timestamp}/{input}"+".tsv"
run_prep_df.to_csv(input,sep ='\t',header=False,index=False)


