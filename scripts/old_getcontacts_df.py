import pandas as pd
import sys
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", required=True,
        help="Input csv file")
parser.add_argument("-v", "--variant_hash", required=False,
        help="Alphanumeric character STRING linked to a specific amino acid sequence")
parser.add_argument("-t", "--timestamp_version", required=False,
        help="System TIMESTAMP applied to distinct analytical applications")
args = parser.parse_args()
print(args.input_file)
input=(args.input_file)
pdb_path = os.path.dirname(args.input_file)

def prepare_(test_file):
    with open(test_file) as parsed_file:
        parsed_file2 = pd.read_csv(parsed_file,low_memory=False,sep="\t",header=None)
        #print(parsed_file)
        print(parsed_file2)
        gc_df = pd.DataFrame(data=parsed_file2).drop(columns=[0])
        #, columns=["index1","interaction_type","atom_1","atom_2"]).astype("str")
        #gc_df.info(verbose=True)
        #print(gc_df)
    return gc_df
run_prep_df = prepare_(input)
variant_hash = args.variant_hash
timestamp_version = args.timestamp_version

if variant_hash:
    run_prep_df[11] = variant_hash
if timestamp_version:
    run_prep_df[12] = pd.to_datetime(timestamp_version)
    
print(run_prep_df)
#structure_df = format_df_for_output(run_prep_df)
filename=pdb_path + "/"+ variant_hash+'.getcontacts.csv'
print(filename)
run_prep_df.to_csv(filename,sep ='\t', header=False, index=False)