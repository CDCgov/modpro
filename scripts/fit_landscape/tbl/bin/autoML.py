#!/usr/bin/env python
# coding: utf-8
# Implementation of autoML Python package to identify "optimal" ML models for applied XYZ+ coordinate factors
# Author: Brian Mann
# Last Modified: 2024-03-26

### IMPORT PYTHON PACKAGES ###

import os
import sys
import re
import argparse
from itertools import repeat

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import plotly.express as px

import pycaret
from pycaret.regression import *


### ASSIGN GLOBAL VARIABLES ###
infile = ''
outfile = ''
factor_context = 'composite'
zcontext = 'all' # Options include: "all", "antigenic", and "proteomic"
cv_fold = 5
folde_iterations = 200


### INPUT VARIABLES ###
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
parser.add_argument(
    "-x",
    "--dependent",
    required=False,
    help="Provide context for applied XY-coordinate axes (i.e., factors) [Options: composite or iterative]"
)
parser.add_argument(
    "-y",
    "--independent",
    required=False,
    help="Provide context for applied Z-coordinate axes (i.e., independent variables) [Options: all, antigenic, or proteomic]"
)

args = parser.parse_args()

infile = args.infile
outfile = args.outfile
# outfile = '/home/ytw4/dev/fit_landscape/tbl/temp/autoML_xyzsurvey_20230326_TestE.txt'
factor_context = args.dependent
zcontext = args.independent


### IMPORT DATA ###
if args.delimiter is None:
    infile = pd.read_csv(infile)
elif args.delimiter == "comma":
    infile = pd.read_csv(infile, sep=",")
elif args.delimiter == "tab":
    infile = pd.read_csv(infile, sep="\t")
elif args.delimiter == "space":
    infile = pd.read_csv(infile, sep=" ")

# infile = pd.read_csv(
#     '/home/ytw4/dev/fit_landscape/tbl/temp/H1_HA_consolidated_impute_antigenicity_VALIDATION_20240228.txt',
#     sep = '\t'
# )


### DEFINITIONS ###
def model_selection( df ):

    # Execute simple autoML regression analysis
    # NOTE: Skip-over all analyses with less than 10 rows with no NaN data entries
    # Why? In order to perform cross-validation, a minimum of [5] folds with [2] data entries is required
    if n_rows > 9:

        # Establish core parameters for 'autoML' package applications with the 'setup()' function
        auto_reg = setup(
            df,
            preprocess=False,
            train_size=0.75,
            imputation_type=None,
            remove_multicollinearity=True,
            remove_outliers=True,
            transformation=True
        )

        # Evaluate applicability of all encoded regression models with a 75/25 "train"/"test" data split using the 'compare_models()' function
        model_selection = auto_reg.compare_models(
            cross_validation=False
        )
                
        # 'pull()' function extracts the output of most recent 'autoML' package application
        results = pull().reset_index()
        iteration_output.append(n_rows)

        # Identify "optimal" ML method based on the TOP-ranked model (i.e., top-listed row)
        optimal_method = results['index'][0]
        iteration_output.append(optimal_method)
        iteration_output.append(results['Model'][0])

        # Extract evaluation metrics for both the TRAIN/TEST models
        optimized_model = create_model(
            optimal_method,
            return_train_score=True,
            cross_validation=False
        )
        results = pull().reset_index()

        # TEST dataset output metrics/analytics
        iteration_output.append(results['MAE'][0])
        iteration_output.append(results['MSE'][0])
        iteration_output.append(results['RMSE'][0])
        iteration_output.append(results['R2'][0])
        iteration_output.append(results['RMSLE'][0])
        iteration_output.append(results['MAPE'][0])

        # TRAIN dataset output metrics/analytics
        iteration_output.append(results['MAE'][1])
        iteration_output.append(results['MSE'][1])
        iteration_output.append(results['RMSE'][1])
        iteration_output.append(results['R2'][1])
        iteration_output.append(results['RMSLE'][1])
        iteration_output.append(results['MAPE'][1])

        # Fill-in NaNs for "skipped" comparisons (due to too few data points)
    else:
        iteration_output.append(n_rows)
        for _ in repeat(None, 14):
            iteration_output.append("NaN")
    
    return()

### MACHINE-LEARNING MODEL SELECTION ###

# Remove uninformative columns from input dataframe
# xcolumns = [
#     "subtype",
#     "protein",
#     "clade",
#     "substitution_list",
#     "n",
#     "surveillance_period",
#     "aa_seq",
#     "data_source",
#     "version"
# ]
# processed = infile.drop(columns=xcolumns)
processed = infile.copy()

# Initiate a lists to store all parameters in each autoML iterations
output_archive = []

# Define list(s) of all valid X/Y (termed X) and Z factors (termed Y; aka columns)
yfactors_rosetta = [
    "total_score",
    "disulfide_potential",
    "attraction_potential",
    "potential_energy",
    "intraresidue_repulsion",
    "weighted_anisotropic_solvation",
    "repulsive_potential",
    "solvation_energy",
    "hbond_bb_sidechain",
    "hbond_long_backbone",
    "hbond_sidechain",
    "hbond_short_backbone",
    "omega_dihedral",
    "tyrosine_torsion"
]

processed_columns = list(processed)
yfactors_antigenic = list(
    filter(
        # lambda x: re.match(r'.*(log).*(cladei).*', x),
        lambda x: re.match(r'.*(log).*', x),
        processed_columns
    )
)

if zcontext == 'all':
    yfactors = yfactors_rosetta + yfactors_antigenic
elif zcontext == 'antigenic':
    yfactors = yfactors_antigenic
elif zcontext == 'proteomic':
    yfactors = yfactors_rosetta

# xfactors = [
#         "hd",
#         "pcd",
#         "carbon_sd_reference",
#         "carbon_sd_seasonal",
#         "calpha_sd_reference",
#         "calpha_sd_seasonal",
#         "heavy_sd_reference",
#         "heavy_sd_seasonal"
#     ]
xfactors = list(
    filter(
        lambda x: re.match(r'.*(hd|pcd|sd).*', x),
        processed_columns
    )
)

# Execute nested FOR loops based on user-defined X factor application(s)
if factor_context == 'composite':
    print(
        "Applied metrics (XY): ",
        xfactors
    )

    # Define Z-coordinate axis from available independent variables
    for k in range(len(yfactors)):
        iteration_output = []
            
        # Print status interval for current Z-coordinate axis
        print(
            "Iteration initiated: ",
            yfactors[k],
            " (Z)"
        )
            
        # Define dataframe with single independent variable (Y factor)
        factor_subset = xfactors.copy()
        factor_subset.append(yfactors[k])
        iteration_output.append(yfactors[k])
    
        # Sub-select processed dataframe to include only selected factors and independent variable in current iteration
        ingest = processed[factor_subset]
    
        # Drop all rows with NaNs
        ingest = ingest.dropna()
    
        # COUNT() number of remaining rows (i.e., variants)
        n_rows = ingest.shape[0]

        model_selection(ingest)
        output_archive.append(iteration_output)

    # Transform autoML results into a dataframe
    output = pd.DataFrame(output_archive)

    # Assign column headers to output dataframe
    output = output.rename(
        columns={
            0: 'z',
            1: 'n',
            2: 'model_abbreviation',
            3: 'model',        
            4: 'mae',
            5: 'mse',
            6: 'rmse',
            7: 'r2',
            8: 'rmsle',
            9: 'mape',
            10: 'mae_train',
            11: 'mse_train',
            12: 'rmse_train',
            13: 'r2_train',
            14: 'rmsle_train',
            15: 'mape_train'
        }
    )

elif factor_context == 'iterative':
    xfactors2 = xfactors.copy()
    xfactors_loo = xfactors.copy()

    # Define X-coordinate axis from available dependent variables
    for i in range(len(xfactors)):
        xfactors_loo.remove(xfactors[i])
    
        # Define Y-coordinate axis from available dependent variables (minus factor assigned in the X-coordinate axis)
        for j in range(len(xfactors_loo)):
        
            # Define Z-coordinate axis from available independent variables
            for k in range(len(yfactors)):
                iteration_output = []
            
                # Print status interval for current XYZ coordinate-axes
                print(
                    "Iteration initiated: ",
                    xfactors[i],
                    " (X), ",
                    xfactors_loo[j],
                    " (Y), and ",
                    yfactors[k],
                    " (Z)"
                )
            
                # Assign both selected X/Y dependent variables as input "X" factors
                xfactors_xy = []
                xfactors_xy.append(xfactors[i])
                xfactors_xy.append(xfactors_loo[j])            
        
                # Define dataframe with all valid dependent variables (X factors) and single independent variable (Y factor)
                factor_subset = xfactors_xy
                factor_subset.append(yfactors[k])
                iteration_output.append(xfactors[i])
                iteration_output.append(xfactors_loo[j])
                iteration_output.append(yfactors[k])
    
                # Sub-select processed dataframe to include only selected factors and independent variable in current iteration
                ingest = processed[factor_subset]
    
                # Drop all rows with NaNs
                ingest = ingest.dropna()
    
                # COUNT() number of remaining rows (i.e., variants)
                n_rows = ingest.shape[0]
    
                model_selection( ingest )
                output_archive.append(iteration_output)
    
        # Iteratively remove elements from X factor list to avoid symmetric analyses
        xfactors2.remove(xfactors[i])
    
        # Terminate loop if all unique X factor pairs have been evaluated
        if len(xfactors2) == 1:
            break

    # Transform autoML results into a dataframe
    output = pd.DataFrame(output_archive)

    # Assign column headers to output dataframe
    output = output.rename(
        columns={
            0: 'x',
            1: 'y',
            2: 'z',
            3: 'n',
            4: 'model_abbreviation',
            5: 'model',        
            6: 'mae',
            7: 'mse',
            8: 'rmse',
            9: 'r2',
            10: 'rmsle',
            11: 'mape',
            12: 'mae_train',
            13: 'mse_train',
            14: 'rmse_train',
            15: 'r2_train',
            16: 'rmsle_train',
            17: 'mape_train'
        }
    )

# Export labeled dataframe to indicated filename/location (as tab-delimited TXT file)
with open(outfile, 'w') as file:
    output.to_csv(
        outfile,
        sep ='\t',
        index=False
    )