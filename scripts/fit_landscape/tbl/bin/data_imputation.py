#!/usr/bin/env python3.8
# Author: Brian Mann (Lead Developer)
# Purpose: Impute (or replace) "missing" data (i.e., NaNs) in provided data matrices
# Created: 2024-07-23
# Input: python3.8 <PYTHON SCRIPT> -i <INPUT FILENAME> -d <DELIMITER> -o <OUTPUT FILENAME>

### IMPORT MODULES ###
import numpy as np
import pandas as pd
import argparse


### DEFINITIONS ###

# System variable definitions/setup: `argparse`
infile = None
outfile = None
delimiter = None

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--infile",
    required=True,
    help="Provide filename for input data matrix"
)
parser.add_argument(
    "-o",
    "--outfile",
    required=True,
    help="Provide filename for output data matrix"
)
parser.add_argument(
    "-d",
    "--delimiter",
    required=False,
    help="Indicate delimiter (and actual character) for provided data file. [Options: comma, tab, and space]"
)

args = parser.parse_args()

infile = args.infile
outfile = args.outfile

if args.delimiter is not None:
    delimiter = args.delimiter


### INGEST DATA ###
if args.delimiter is None:
    ingest = pd.read_csv(infile)
elif args.delimiter == "comma":
    ingest = pd.read_csv(infile, sep=",")
elif args.delimiter == "tab":
    ingest = pd.read_csv(infile, sep="\t")
elif args.delimiter == "space":
    ingest = pd.read_csv(infile, sep=" ")

# Extract all columns with antigenic data in addition to the "clade" and "variant_hash" identifier columns
ingest = ingest.replace('*',np.nan)
ingest = ingest.drop_duplicates()
antigenic = ingest.filter(regex='clade|log|titer|ag\_cdc\_id|ag\_isolate\_id|variant\_hash')

# Drop all duplicate rows from antigenically-subset dataframe (Rationale: All 'variant_hash' entries should be unique)
antigenic = antigenic.drop_duplicates().reset_index()
antigenic = antigenic.filter(regex='clade|log|titer|ag\_cdc\_id|ag\_isolate\_id|variant\_hash')
antigenic_parse = antigenic.copy()


### PROCESS DATA ###
def log_transform( titer ):
    if isinstance( titer, (int, float)):
        titer_adjusted = titer / 5
        i = np.log2( titer_adjusted )
    else:
        i = np.nan
    return i

def back_transform( titerlog ):
    if isinstance( titerlog, (int, float)):
        i = pow( 2, titerlog ) * 5
    else:
        i = np.nan
    return i

# Re-order "clade" and "variant_hash" as leading columns (required for column indexing in following sub-routines)
antigenic_parse = antigenic_parse.drop(columns=['clade'])
antigenic_parse.insert(loc=0, column='clade', value=antigenic['clade'])
antigenic_parse = antigenic_parse.drop(columns=['variant_hash'])
antigenic_parse.insert(loc=0, column='variant_hash', value=antigenic['variant_hash'])

# Format all antigenic data columns as the FLOAT data type
antigenic_parse = antigenic_parse.astype({col: float for col in antigenic_parse.columns[2:]})

# Normalize all HI and/or HINT raw titer values in LOG(2)-space
cols = list(filter(lambda i: 'titer' in i, antigenic_parse.columns))
for i in cols:
    antigenic_parse[i] = antigenic_parse[i].apply(log_transform)

# Compute clade-specific, mean log-transformed titer values (per column)
antigenic_clade = antigenic_parse.filter(regex='clade|log|titer')
clade_mean_gmt = antigenic_clade.groupby('clade').mean(numeric_only=True).reset_index()

# Replace all NAs in dataframe with geometric mean titer (GMT) per clade [NOTE: NAs may persist based on data availability]
antigenic_cladei = antigenic_parse.copy()
cols = list(filter(lambda i: 'logfold' in i, antigenic_cladei.columns))
cols_titer = list(filter(lambda i: 'titer' in i, antigenic_cladei.columns))
for i in cols:
    j = i + "_y"
    k = i + "|cladei"
    antigenic_cladei[i] = antigenic_cladei[i].fillna(antigenic_cladei.merge(clade_mean_gmt, on='clade', how='left')[j])
    if i in cols_titer:
        antigenic_cladei[i] = antigenic_cladei[i].apply(back_transform)
    antigenic_cladei = antigenic_cladei.rename(columns={i: k})

# Replace all NAs in dataframe with GMT of the overall column
antigenic_meani = antigenic_parse.copy()
cols = list(filter(lambda i: 'logfold' in i, antigenic_meani.columns))
cols_titer = list(filter(lambda i: 'titer' in i, antigenic_meani.columns))
for i in cols:
    j = i + "|meani"
    antigenic_meani[i] = antigenic_meani[i].fillna(antigenic_meani.mean(numeric_only=True)[i])
    if i in cols_titer:
        antigenic_meani[i] = antigenic_meani[i].apply(back_transform)
    antigenic_meani = antigenic_meani.rename(columns={i: j})

# Merge all imputed data into original dataframe [NOTE: New columns annotated with additional suffix(es)]
cols = list(ingest.filter(regex='variant\_hash|clade|ag\_cdc\_id|ag\_isolate\_id').columns)
consolidated = ingest.merge(antigenic_cladei, on=cols, how='left')
consolidated = consolidated.merge(antigenic_meani, on=cols, how='left')

# Write results to user-defined file in tab-delimited format (if file already exists)
with open(outfile, 'w') as file:
    consolidated.to_csv(outfile, sep ='\t', index=False)