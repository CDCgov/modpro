#!/usr/bin/env python
# coding: utf-8
# Purpose: Compute number of "structural differences" between user-defined and "reference" PDB protein structures
# Concept URL: https://pubmed.ncbi.nlm.nih.gov/23277561/
# Author: Brian Mann (Lead Developer)
# Last Modified: 2025-05-07

### IMPORT PYTHON PACKAGES ###

import os
import sys
import re
import argparse
import scipy as sp
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

### ASSIGN GLOBAL VARIABLES ###

# System variable definitions/setup: `argparse`
Input = None
Output = None
Context = None
Relaxed = False
Limit = 4.5
Confidence = 70.0
Proteomic_context = None
Site_range = None


### INPUT VARIABLES ###

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--input",
    required=True,
    help="Provide input PDB-formatted file"
)
parser.add_argument(
    "-o",
    "--output",
    required=True,
    help="Provide filename and location for output tab-delimited TXT file"
)
parser.add_argument(
    "-c",
    "--context",
    required=False,
    help="Indicate context for structural distance computation: 'calpha', 'carbon', or 'heavy' [Default: 'calpha']"
)
parser.add_argument(
    "-l",
    "--limit",
    required=False,
    help="Indicate maximal distance between atoms (in Angstroms) [Default: 4.5]"
)
parser.add_argument(
    "-k",
    "--confidence",
    required=False,
    help="Indicate minimum residue pLDDT value [Default: 70.0]"
)
parser.add_argument(
    "-x",
    "--proteomic_context",
    required=False,
    help="Indicate context for user-defined amino acid positions in applied computations [Default: HA; Examples: HA, HA1, functional_sites, etc.]"
)
parser.add_argument(
    "-r",
    "--relaxed",
    required=False,
    help="Provide indication True or False whether the applied PDB file has been 'relaxed' by the Rosetta software package"
)
parser.add_argument(
    "-s",
    "--site_range",
    required=False,
    help="Indicate semicolon-delimited list of valid amino acid positions; annotate ranges with '..' characters"
)

args = parser.parse_args()

Input = args.input
Output = args.output
valid_relaxed_entries = [
    True,
    'T',
    'true',
    'True'
]

if args.context is not None:
    Context = args.context
if args.limit is not None:
    Limit = float(args.limit)
if args.confidence is not None:
    Confidence = float(args.confidence)
if args.proteomic_context is not None:
    Proteomic_context = args.proteomic_context
if args.relaxed is not None:
    if args.relaxed in valid_relaxed_entries:
        Relaxed = True
        Confidence = None
if args.site_range is not None:
    Site_range = args.site_range


### DEFINE FUNCTIONS ###

# Parse ingested PDB files into dataframes with X-, Y-, and Z-coordinates defined
# Remove all non-"ATOM" rows (retain only atomic coordinates)
# Add extra leading space (\s) to all "-" characters - fixes concatenated negative X, Y, and/or Z-coordinates
# Remove trailing whitespace per line (avoids potential NULL trailing columns)
# Extract user-defined AA_ID range(s) (apply the same PDB file to compute "strucutral distance" for HA, HA1, HA2, antibody-binding site, etc.)
def parse_pdb( x ):
    pdb_parsed = x[~x[0].str.contains('|'.join(["TER","END"]))] \
        .replace('-',' -',regex=True) \
        .replace(r" +$", r"",regex=True)
    pdb_split = pdb_parsed[0].str.split('\s+',expand=True)
    pdb_split.columns = [
        'index_type',
        'atom_id',
        'element_id',
        'aa',
        'chain_id',
        'aa_id',
        'x',
        'y',
        'z',
        'i',
        'plddt',
        'element'
    ]
    pdb_split = pdb_split.astype(
        dtype= {
            "atom_id":"int",
            "aa_id":"int",
            "x":"float64",
            "y":"float64",
            "z":"float64",
            "plddt":"float64"
        }
    )
    pdb_split = pdb_split[pdb_split['element_id'].isin(atomic_target)]
    if len(sfilter) > 0:
        pdb_split = pdb_split[pdb_split['aa_id'].isin(sfilter)]
    return(pdb_split)

def structural_contacts(x):
    # Exclude all {XYZ} coordinates with "low" confidence: pLDDT < user-defined 'confidence' variable
    # Replace all "low" confidence {XYZ} coordinates with NaN entries
    for i in range(len(coordinates)):
        if Relaxed == False:
            x.loc[x['plddt']<Confidence,coordinates[i]] = np.nan
    # Extract {XYZ} coordinate columns to compute distance matrix
    # URL: https://www.dabblingbadger.com/blog/2020/2/27/implementing-euclidean-distance-matrix-calculations-from-scratch-in-python
    x_array = x[coordinates].reset_index(drop=True).to_numpy()
    # Compute full "intra"-distance matrix of all retained atomic {XYZ} coordinates for both the query (q) and reference (r) PDB files
    xd_array = sp.spatial.distance_matrix(x_array,x_array)
    # Replace all array entries with distance >= user-provided bond "limit" (e.g., 4.5) to "0"
    # Replace all remainder array entries with distance > 0 to "1"
    # Render all NaN entries as "0" via Numpy `nan_to_num()` function
    xd_array[xd_array>=Limit] = 0
    xd_array[xd_array>0] = 1
    return(xd_array)

def adjust_index(x,y,z):
    # Standardize array shape with NULL entries to cover discrepant amino acid sequence lengths (i.e., cross-subtype)
    i = list(np.arange(((x*z)-1),((x*z)-(1/(y+1))),(1/y)))
    i = i + list(np.arange(((x*z*2)+y-1),((x*z*2)+y-(1/(y+1))),(1/y)))
    i = i + list(np.arange(((x*z*3)+(y*2)-1),((x*z*3)+(y*2)-(1/(y+1))),(1/y)))
    i = np.add(i,(1/(y+2)))
    return(i)


### PARAMETER/VARIABLE VALIDATION ###

# Extract identifier from user-provided filename [Example: variant hash]
# Detection of Rosetta "relaxed" PDB structures based on the presence/absence of the "relaxed" or "relaxed_quick_" input file prefix(es)
if '.pdb' not in str(Input):
    print(
        Input,
        "\nInvalid argument for --input; -i parameter\nPDB file format (.pdb) required"
    )
    sys.exit(2)
else:
    Input_id = re.search('(^.*)\.pdb',re.split('/',Input)[-1]).group(1)
    if "relax" in Input_id:
        Relaxed = True
        Input_id = re.split('_',re.split('relaxed_quick_|relaxed_',Input_id)[1])[0]

# Which elements in the applied PDB file will be applied in the Euclidean distance computaions?
if Context == "calpha":
    atomic_target = ["CA"]
elif Context == "carbon":
    atomic_target = ["CA","C"]
elif Context == "heavy":
    atomic_target = ["CA","C","N","O"]
else:
    print(
        "\nInvalid argument for --context; -c parameter\nOptions: calpha, carbon, or heavy"
    )
    sys.exit(2)

# Verify user-defined structural distance "limit" is in INT or FLOAT format
if isinstance(Limit,float) | isinstance(Limit,int):
    Limit = float(Limit)
else:
    print(
        "\nInvalid argument for --limit; -l parameter\nInteger or decimal value required"
    )
    sys.exit(2)

# Verify user-defined limit atomic positional "confidence" (i.e., pLDDT value in AlphaFold2) is in INT or FLOAT format
# NOTE: "Relaxed" PDB structures optimized by the Rosetta software suite will nullify applicability of the pLDDT column with all values set to 0.0
# Implementation of "relaxed" PDB structures requires inclusion of the additional --confidence input argument with an entry of "0.0" 
if isinstance(Confidence,float) | isinstance(Confidence,int):
    Confidence = float(Confidence)
else:
    if Confidence != None:
        print(
            "\nInvalid argument for --confidence; -k parameter\nInteger or decimal value required"
        )
        sys.exit(2)

# Expand user-provided ";" and ".." delimited list of applicable amino acid positions
# URL: https://stackoverflow.com/questions/30201119/expanding-a-block-of-numbers-in-python
# Example: '95;100..110;120..125;135' [site_range input]
sfilter = []
if Site_range is not None:
    srange_list = list(Site_range.split(';'))
    
    for i in srange_list:
        if '..' in i:
            Start, Stop = map(int,i.split('..'))
            sfilter.extend(range(Start,(Stop+1)))
        else:
            sfilter.append(int(i))

sfilter.sort()


# Define all "reference" X-, Y-, and Z-coordinates (per subtype/context)
# NOTE: "Relaxed" reference PDB structures are applied if relaxed PDB files are processed
if Relaxed:
    pdb_reference = [
        'relaxed_B_HA_reference',
        'relaxed_B_HA_seasonal',
        'relaxed_H1_HA_reference',
        'relaxed_H1_HA_seasonal',
        'relaxed_H3_HA_reference',
        'relaxed_H3_HA_seasonal'
    ]
else:
    pdb_reference = [
        'B_HA_reference',
        'B_HA_seasonal',
        'H1_HA_reference',
        'H1_HA_seasonal',
        'H3_HA_reference',
        'H3_HA_seasonal'
    ]

# Define valid coordinate axes in Euclidean distance matrices
coordinates = [
    'x',
    'y',
    'z'
]


### EXECUTE CODE ###

# Parse user-define input PDB file into atomic-level dataframe
pdb = pd.read_csv( 
    Input, 
    header=None,
    sep='\t'
)

# Relaxed PDB structures processed by the Rosetta software package include additional header/footer metadata which need to be removed to ensure a standardized dataframe shape
if Relaxed:
    pdb = pdb[pdb[0].str.contains('ATOM|TER')].reset_index(drop=True)
    

# Implement 'for' loop over all defined historical and seasonal "references" per subtype to compute structural distance (SD) metric
sd_results = []

for reference in pdb_reference:
    query_pdb = parse_pdb(pdb).reset_index(drop=True)
    
    Ref_input = ROOT / "scripts"/ "tbl" / f"{reference}.pdb"
    Ref_label = f"{reference}_pdb"
    Labels = reference.split("_")
    ref_pdb = pd.read_csv(
        Ref_input,
        header=None,
        sep='\t'
    )
    if Relaxed:
        del Labels[0]
        ref_pdb = ref_pdb[ref_pdb[0].str.contains('ATOM|TER')].reset_index(drop=True)
    reference_pdb = parse_pdb(ref_pdb).reset_index(drop=True)
    
    if reference_pdb.shape[0] != query_pdb.shape[0]:
        Discrepant_rows = abs((query_pdb.shape[0]-reference_pdb.shape[0])/(reference_pdb['chain_id'].nunique()))#*len(atomic_target)))
        row_insertion = np.empty([int(Discrepant_rows),int(reference_pdb.shape[1])])
        row_insertion[:] = np.nan

        if reference_pdb.shape[0] > query_pdb.shape[0]:
            index_iter = max(query_pdb['aa_id'])
            index_insert = adjust_index(index_iter,Discrepant_rows,len(atomic_target))
            for i in index_insert:
                query_pdb.loc[i] = row_insertion[0]
                query_pdb = query_pdb.sort_index().reset_index(drop=True)
        else:
            index_iter = max(reference_pdb['aa_id'])
            index_insert = adjust_index(index_iter,Discrepant_rows,len(atomic_target))
            for i in index_insert:
                reference_pdb.loc[i] = row_insertion[0]
                reference_pdb = reference_pdb.sort_index().reset_index(drop=True)
    
    # Identify all molecular bonds in applied PDB file with a bond-length <= user-defined limit [Default: 0.45nm, 4.5 Angstroms]
    qd_array = structural_contacts(query_pdb)
    rd_array = structural_contacts(reference_pdb)
    
    if rd_array.shape == qd_array.shape:
        diff_array = abs(rd_array-qd_array)
        Sd_metric = int(np.sum(np.nan_to_num(diff_array))/2)
        sd_list = [
            Labels[1],
            Labels[0],
            Labels[2],
            Input_id,
            Sd_metric
        ]
        sd_results.append(sd_list)
    else:
        continue


### WRITE-TO-FILE ###

# Process results into exportable dataframe/data-types
outfile = pd.DataFrame(sd_results)
outfile.insert(
    4,
    "sd_context",
    Context,
    True
)
outfile.insert(
    5,
    "sd_limit",
    Limit,
    True
)
outfile.insert(
    6,
    "sd_confidence",
    Confidence,
    True
)
outfile.insert(
    7,
    "proteomic_context",
    Proteomic_context,
    True
)
outfile.insert(
    8,
    "site_range",
    Site_range,
    True
)
outfile.columns = [
    'protein',
    'subtype',
    'reference_context',
    'identifier',
    'sd_context',
    'sd_limit',
    'sd_confidence',
    'proteomic_context',
    'site_range',
    'sd'
]
outfile = outfile.astype(
    dtype= {
        "protein":"string",
        "subtype":"string",
        "reference_context":"string",
        "identifier":"string",
        "sd_context":"string",
        "sd_limit":"float64",
        "sd_confidence":"float64",
        "proteomic_context":"string",
        "site_range":"string",
        "sd":"int"
    }
)

# Write (or append) results to user-defined file in tab-delimited format (if file already exists)
#if os.path.exists(Output):
    #with open(Output,'a') as file:
        #outfile.to_csv(
            #file,
            #sep ='\t',
            #header=False,
            #index=False
        #)
#else:
with open(Output,'w') as file:
    outfile.to_csv(
        file,
        sep ='\t',
        header=False,
        index=False
    )

