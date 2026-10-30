#!/usr/bin/env python
# coding: utf-8
# Merge all data files as staging for autoML and maching-learning analytics
# Author: Brian Mann
# Last Modified: 2024-07-22

### IMPORT PYTHON PACKAGES ###

import re
import argparse
import pandas as pd


### ASSIGN GLOBAL VARIABLES ###

antigenic = None
distance = None
proteomic = None
individual = False


### INPUT VARIABLES ###

parser = argparse.ArgumentParser()

parser.add_argument(
    "-a",
    "--antigenic",
    required=False,
    help="Provide input file containing antigenic trends."
)

parser.add_argument(
    "-d",
    "--distance",
    required=False,
    help="Provide input file containing structural distance computations."
)

parser.add_argument(
    "-p",
    "--proteomic",
    required=False,
    help="Provide input file containing Rosetta energy computations."
)

parser.add_argument(
    "-i",
    "--individual",
    required=False,
    help="Provide BOOLEAN (True) confirmation if antigenic data is aggregated by individual antigen ISOLATE IDs."
)

args = parser.parse_args()

if args.antigenic is not None:
    antigenic = args.antigenic
if args.distance is not None:
    distance = args.distance
if args.proteomic is not None:
    proteomic = args.proteomic
if args.individual is not None:
    individual = args.individual


### IMPORT DATA ###

if args.antigenic is None and individual:
    antigenic_df = pd.read_csv('tbl/indv_antigenic_summaries.txt',sep='\t')
elif args.antigenic is None:
    antigenic_df = pd.read_csv('tbl/antigenic_summaries.txt',sep='\t')
else:
    antigenic_df = pd.read_csv(antigenic)

if args.distance is None:
    distance_df = pd.read_csv('tbl/genetic_summaries.txt',sep='\t')
else:
    distance_df = pd.read_csv(distance)

distance_label_df = distance_df.rename(
    columns={
        0: 'subtype',
        1: 'protein',
        2: 'clade',
        3: 'identifier',
        4: 'hd',
        5: 'pcd',
        6: 'substitution_list',
        7: 'aa_seq'
    }
)

if args.proteomic is None:
    proteomic_df = pd.read_csv( 'tbl/rosetta_summaries.txt', sep = '\t' )
else:
    proteomic_df = pd.read_csv(proteomic)

rosetta_references = pd.read_csv( 'tbl/rosetta_references.txt', sep = '\t' )
calpha_sd_df = pd.read_csv( 'tbl/structural_distance_calpha.txt', sep = '\t' )
carbon_sd_df = pd.read_csv( 'tbl/structural_distance_carbon.txt', sep = '\t' )
heavy_sd_df = pd.read_csv( 'tbl/structural_distance_heavy.txt', sep = '\t' )

proteomic_label_df = proteomic_df.rename(
    columns={
        0: 'total_score',
        1: 'disulfide_potential',
        2: 'attraction_potential',
        3: 'potential_energy',        
        4: 'intraresidue_repulsion',
        5: 'weighted_anisotropic_solvation',
        6: 'repulsive_potential',
        7: 'solvation_energy',
        8: 'hbond_bb_sidechain',
        9: 'hbond_long_bb',
        10: 'hbond_sidechain',
        11: 'hbond_short_bb',
        12: 'omega_dihedral',
        13: 'tyrosine_torsion',
        14: 'identifier',
        15: 'version'
    }
)


### PROCESS INPUT DATA ###

# Designate consolidated column headers for provided antigenic trends
antigenic_df['subtype'] = antigenic_df['subtype'].replace( 'H1 swl', 'H1' )
antigenic_df['subtype'] = antigenic_df['subtype'].replace( 'B vic', 'B' )
antigenic_df['aggregate_label'] = antigenic_df['subtype'].map(str) + '|' + antigenic_df['test_protocol'].map(str) + '|' + antigenic_df['aggregation'].map(str) + '|' + antigenic_df['lot'].map(str) + '|logfold'

# PIVOT antigenic data to summarize trends column-wise (per 'test_protocol'/'lot')
if individual:
    antigenic_subset_df = antigenic_df[[
        'aggregate_label',
        'identifier',
        'ag_cdc_id',
        'ag_isolate_id',
        'avg_logfold'
    ]].sort_values(by=['identifier','ag_cdc_id','ag_isolate_id','aggregate_label']).reset_index(drop=True)
        
    antigenic_pivot_df = antigenic_subset_df.pivot_table(
        index=['identifier','ag_cdc_id','ag_isolate_id'],
        columns='aggregate_label',
        values='avg_logfold'
    ).reset_index()
else:
    antigenic_subset_df = antigenic_df[[
        'aggregate_label',
        'identifier',
        'avg_logfold'
    ]].sort_values(by=['identifier','aggregate_label']).reset_index(drop=True)
    
    antigenic_pivot_df = antigenic_subset_df.pivot_table(
        index='identifier',
        columns='aggregate_label',
        values='avg_logfold'
    ).reset_index()

# PIVOT structural distance to summarize trends column-wise (per 'variant_hash')
contexts = ['calpha','carbon','heavy']
for i in range(len(contexts)):
    if contexts[i] == 'calpha':
        sd_label_df = calpha_sd_df.copy()
    elif contexts[i] == 'carbon':
        sd_label_df = carbon_sd_df.copy()
    elif contexts[i] == 'heavy':
        sd_label_df = heavy_sd_df.copy()
    
    sd_label_df['label'] = sd_label_df['subtype'].map(str) + '_' + sd_label_df['reference_context'].map(str) + '_' + sd_label_df['atomic_context'] + '_sd'
    
    sd_subset_df = sd_label_df[[
        'identifier',
        'label',
        'sd'
    ]].sort_values(by=['identifier','label']).reset_index(drop=True)
    
    sd_pivot_df = sd_subset_df.pivot_table(
        index='identifier',
        columns='label',
        values='sd'
    ).reset_index()
    
    if contexts[i] == 'calpha':
        calpha_sd_pivot_df = sd_pivot_df.copy()
    elif contexts[i] == 'carbon':
        carbon_sd_pivot_df = sd_pivot_df.copy()
    elif contexts[i] == 'heavy':
        heavy_sd_pivot_df = sd_pivot_df.copy()

# Merge all structural distance dataframes (by 'identifier' column)
cc_sd_merge_df = calpha_sd_pivot_df.merge( carbon_sd_pivot_df, on = ['identifier'], how = 'outer' )
sd_merge_df = cc_sd_merge_df.merge( heavy_sd_pivot_df, on = ['identifier'], how = 'outer' )

# Extract "relative" Rosetta metrics based on seasonal/historical references
prosetta_df = rosetta_references.merge( proteomic_label_df, on = ['identifier'], how = 'left' )

prosetta_reference_df = prosetta_df[ prosetta_df['context'] == 'reference' ]
prosetta_seasonal_df = prosetta_df[ prosetta_df['context'] == 'seasonal' ]
prosetta_merge_df = prosetta_reference_df.merge( prosetta_seasonal_df, on = 'subtype', how = 'outer' )

# Merge structural distance and proteomic dataframes
psd_merge_df = sd_merge_df.merge( proteomic_label_df, on = ['identifier'], how = 'outer' )

# Merge combined dataframe with genetic summary
dpsd_merge_df = distance_label_df.merge( psd_merge_df, on = ['identifier'], how = 'outer' )

# Merge "relative" Rosetta dataframes & compute "relative" Rosetta metrics
dpsd_merge_relative_df = dpsd_merge_df.merge( prosetta_merge_df, on = ['subtype'], how = 'left' )

dpsd_merge_relative_df['rrel_total_score'] = dpsd_merge_relative_df['total_score_x'] - dpsd_merge_relative_df['total_score']
dpsd_merge_relative_df['srel_total_score'] = dpsd_merge_relative_df['total_score_y'] - dpsd_merge_relative_df['total_score']

# Merge antigenic and strucutral distance dataframes
merged_df = dpsd_merge_relative_df.merge( antigenic_pivot_df, on = ['identifier'], how = 'outer' )

# Separate merged dataframe into subtype-specific subsets
base_columns = list( distance_label_df )
proteomic_columns = list( proteomic_label_df )
proteomic_subset_columns = proteomic_columns[0:14] + ['rrel_total_score','srel_total_score']
processed_columns = list( sd_merge_df ) + list( antigenic_pivot_df )

if individual:
    specimen_columns = ['ag_cdc_id','ag_isolate_id']

subtypes = ['B vic','H1 swl','H3']
for i in range(len(subtypes)):
    j = subtypes[i]
    if subtypes[i] == 'B vic':
        k = 'B'
    elif subtypes[i] == 'H1 swl':
        k = 'H1'
    else:
        k = subtypes[i]
    
    xtra_columns = list(
        filter(
            lambda x: re.match(r'.*(%s(\_|\|)).*' % k, x),
            processed_columns
        )
    )
    
    if individual:
        columns_subset = base_columns + specimen_columns + proteomic_subset_columns + xtra_columns
        subset_df = merged_df[columns_subset]
    else:
        columns_subset = base_columns + proteomic_subset_columns + xtra_columns
        subset_df = merged_df[columns_subset]

    if subtypes[i] == 'B vic':
        Bvic_df = subset_df[ subset_df['subtype'] == j ].reset_index(drop=True)
    elif subtypes[i] == 'H1 swl':
        H1swl_df = subset_df[ subset_df['subtype'] == j ].reset_index(drop=True)
    elif subtypes[i] == 'H3':
        H3_df = subset_df[ subset_df['subtype'] == j ].reset_index(drop=True)

Bvic_df.rename( { 'identifier' : 'variant_hash' }, axis = 'columns', inplace = True )
H1swl_df.rename( { 'identifier' : 'variant_hash' }, axis = 'columns', inplace = True )
H3_df.rename( { 'identifier' : 'variant_hash' }, axis = 'columns', inplace = True )

### EXPORT PROCESSED DATA ###

if individual:
    Bvic_outfile = 'tbl/indv_B_consolidated.txt'
else:
    Bvic_outfile = 'tbl/B_consolidated.txt'
with open( Bvic_outfile, 'w' ) as file:
    Bvic_df.to_csv(
        Bvic_outfile,
        sep ='\t',
        index = False
    )

if individual:
    H1swl_outfile = 'tbl/indv_H1_consolidated.txt'
else:
    H1swl_outfile = 'tbl/H1_consolidated.txt'
with open( H1swl_outfile, 'w' ) as file:
    H1swl_df.to_csv(
        H1swl_outfile,
        sep ='\t',
        index = False
    )

if individual:
    H3_outfile = 'tbl/indv_H3_consolidated.txt'
else:
    H3_outfile = 'tbl/H3_consolidated.txt'
with open( H3_outfile, 'w' ) as file:
    H3_df.to_csv(
        H3_outfile,
        sep ='\t',
        index = False
    )