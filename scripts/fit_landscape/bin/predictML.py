#!/usr/bin/env python
# coding: utf-8
# Implementation of ensemble ML methods in Python
# Author: Brian Mann
# Last Modified: 2024-07-23

### IMPORT PYTHON PACKAGES ###

import os
import re
import argparse

import numpy as np
import pandas as pd

from sklearn.model_selection import cross_validate
from sklearn.model_selection import RepeatedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import AdaBoostRegressor, GradientBoostingRegressor, ExtraTreesRegressor

### ASSIGN GLOBAL VARIABLES ###

infile = None
outfile = None
delimiter = "tab"
factor_context = "iterative"
zcontext = "mixed" # Options include: "all", "antigenic", "proteomic", and "mixed"
zfilter = None


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
    "-z",
    "--independent",
    required=False,
    help="Provide context for applied Z-coordinate axes (i.e., independent variables) [Options: all, antigenic, or proteomic]"
)

parser.add_argument(
    "-b",
    "--zfilter",
    required=False,
    help="Provide file to filter applied Z-axis coordinates"
)

args = parser.parse_args()

infile = args.infile
outfile = args.outfile
delimiter = args.delimiter
factor_context = args.dependent
zcontext = args.independent
zfilter = args.zfilter


### IMPORT DATA ###

if args.delimiter is None:
    infile = pd.read_csv(infile)
elif args.delimiter == "comma":
    infile = pd.read_csv(infile, sep=",")
elif args.delimiter == "tab":
    infile = pd.read_csv(infile, sep="\t")
elif args.delimiter == "space":
    infile = pd.read_csv(infile, sep=" ")

if args.zfilter is not None:
    zfilter = pd.read_csv(
        args.zfilter,
        sep = '\t'
    )
    seasonal_references = list( zfilter['lot'] )


### DEFINTIONS ###

def dynamic_interval( cmin, cmax ):
    crangei = abs( cmax - cmin ) / 10
    add_rangei = round( cmax + crangei, -int( np.floor( np.log10( crangei ))))
    sub_rangei = round( cmin - crangei, -int( np.floor( np.log10( crangei ))))
    if cmin > 0 and sub_rangei <= 0:
        lbound = 0
        ubound = add_rangei
    elif cmax < 0 and cmax + crangei > 0 :
        lbound = sub_rangei
        ubound = 0
    else:
        lbound = sub_rangei
        ubound = add_rangei
    interval = ( ubound - lbound ) / 25
    meshgrid_axis = np.arange( lbound, ubound, interval )
    return( meshgrid_axis )


### DEFINE X and Y DATAFRAMES ###

processed = infile.copy()
processed_columns = list(processed)
nlimit = int( np.ceil( processed.shape[0] * 0.25 ))

# Define list(s) of all valid X/Y and Z factors
zfactors_rosetta = [
    "total_score",
    "rrel_total_score",
    "srel_total_score",
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

zfactors_antigenic = list(
    filter(
        lambda x: re.match(r'.*(logfold)($|\|cladei$)', x),
        processed_columns
    )
)

if args.zfilter is not None:
    zfactors_antigenic = list(
        filter(
            lambda x: any( substring in x for substring in seasonal_references ),
            zfactors_antigenic
        )
    )

if zcontext == 'all':
    zfactors = zfactors_rosetta + zfactors_antigenic
elif zcontext == 'antigenic':
    zfactors = zfactors_antigenic
elif zcontext == 'mixed':
    #zfactors = ['total_score','rrel_total_score','srel_total_score'] + zfactors_antigenic
    zfactors = ['srel_total_score'] + zfactors_antigenic
elif zcontext == 'proteomic':
    zfactors = zfactors_rosetta

xfactors = list(
    filter(
        lambda x: re.match(r'.*(hd|pcd).*', x),
        processed_columns
    )
)

yfactors = list(
    filter(
        # lambda x: re.match(r'.*(sd|total\_score).*', x),
        lambda x: re.match(r'.*(srel|seasonal).*(sd|total\_score).*', x),
        processed_columns
    )
)

metadata = list(
    filter(
        lambda x: re.match(r'(variant|clade|subtype|protein|substitution).*', x),
        processed_columns
    )
)

pred_landscape_template = pd.DataFrame(
    columns = [
        'subtype',
        'protein',
        'x_factor',
        'y_factor',
        'z_factor',
        'model',
        'n_estimators',
        'learning_rate',
        'max_features',
        'min_samples_split',
        'mae_test',
        'mae_train',
        'r2_test',
        'r2_train',
        'cv_mae_mean',
        'cv_mae_median',
        'cv_r2_mean',
        'cv_r2_median',
        'x_predicted',
        'y_predicted',
        'z_predicted'
    ]
)

summary_outfile = outfile.replace( '.txt', '_summary.txt' )
pred_outfile = outfile.replace( '.txt', '_metrics.txt')

for i in range(len(xfactors)):
    for j in range(len(yfactors)):
        for k in range(len(zfactors)):
            
            pred_landscape_iteration = pred_landscape_template.copy()
            xyz_factors = []
            xy_factors = []
            composite_factors = metadata[0:5]
            
            xyz_factors.append(xfactors[i])
            xyz_factors.append(yfactors[j])
            xyz_factors.append(zfactors[k])
            if len( pd.unique( xyz_factors )) < 3:
                continue
            
            pcomposite_factors = composite_factors.copy()
            pcomposite_factors.extend([xfactors[i],yfactors[j]])

            composite_factors.extend(xyz_factors[0:3])
            xy_factors.append(xfactors[i])
            xy_factors.append(yfactors[j])

            XYZ = processed[xyz_factors].dropna().reset_index(drop=True)
            composite = processed[composite_factors].dropna().reset_index(drop=True)
            pcomposite = processed[pcomposite_factors].dropna().drop_duplicates().reset_index(drop=True)
            pXY = pcomposite[xy_factors].dropna().drop_duplicates().reset_index(drop=True)
            
            if composite.shape[0] < nlimit:
                continue
            
            print( 
                'Iteration Initiated:',
                xfactors[i],
                '(X-axis);',
                yfactors[j],
                '(Y-axis);',
                zfactors[k],
                '(Z-axis)'
            )
            
            X = XYZ[xy_factors]
            Y = XYZ[zfactors[k]]
            
            imin = np.min(X[xfactors[i]])
            imax = np.max(X[xfactors[i]])
            jmin = np.min(X[yfactors[j]])
            jmax = np.max(X[yfactors[j]])
            
            # Interleave two arrays into a 2D matrix with self-transposition of "meshed" 2D matrix in an XY-coordinate format
            # https://numpy.org/doc/stable/reference/generated/numpy.ndarray.T.html (inspired by Bing CoPilot search)
            X_meshgrid = np.meshgrid(
                dynamic_interval( imin, imax ),
                dynamic_interval( jmin, jmax )
            )
            
            X_Predict = pd.DataFrame(
                np.array( X_meshgrid ).T.reshape(-1, 2)
            )
            
            ### RUN MACHINE-LEARNING REGRESSOR ###
            models = [
                AdaBoostRegressor(),
                GradientBoostingRegressor(),
                ExtraTreesRegressor()
            ]
            
            cv_archive = []
            for model in models:
                cv_iteration = []
                print( "Initiate", model, "model ...")
                
                # Scan grid-space to optimize model parameters/hyperparameters
                # Assign model parameter/hyperparameter ranges or limits 
                grid = dict()
                if str(model) in ('AdaBoostRegressor()','GradientBoostingRegressor()'):
                    grid['n_estimators'] = [10, 50, 100, 200]
                    grid['learning_rate'] = [0.0001, 0.001, 0.01, 0.1, 1.0]
                elif str(model) in ('ExtraTreesRegressor()'):
                    grid['n_estimators'] = [10, 50, 100, 200]
                    grid['max_features'] = [1, 2, 3, 4, 5, 6, 7, 8]
                    grid['min_samples_split'] = [2, 3, 4, 5]
                
                cross_validation = RepeatedKFold(
                    n_splits = 5,
                    n_repeats = 20,
                    random_state = 12345678
                )
                
                # Multiple sources cited to observe step-wise grid search optimization
                # URL: https://stackoverflow.com/questions/53973563/using-multiple-metric-evaluation-with-gridsearchcv (view multiple metrics)
                grid_search = GridSearchCV(
                    estimator = model,
                    param_grid = grid,
                    n_jobs = -1,
                    cv = cross_validation,
                    scoring = [
                        'r2',
                        'neg_mean_absolute_error'
                    ],
                    refit = 'neg_mean_absolute_error',
                    return_train_score = True,
                    verbose = 1
                )
                
                grid_result = grid_search.fit(
                    X,
                    Y
                )
                
                if str(model) in ('AdaBoostRegressor()','GradientBoostingRegressor()'):
                    n_estimators = int(grid_result.best_params_['n_estimators'])
                    learning_rate = float(grid_result.best_params_['learning_rate'])
                    max_features = None
                    min_samples_split = None
                    print( 
                        "Optimal Hyperparameters:",
                        n_estimators, "(Estimators) and",
                        learning_rate, "(Learning Rate)"
                    )
                elif str(model) in ('ExtraTreesRegressor()'):
                    n_estimators = int(grid_result.best_params_['n_estimators'])
                    learning_rate = None
                    max_features = int(grid_result.best_params_['max_features'])
                    min_samples_split = int(grid_result.best_params_['min_samples_split'])
                    print( 
                        "Optimal Hyperparameters:",
                        n_estimators, "(Estimators),",
                        max_features, "(Max. Features), and",
                        min_samples_split, "(Min. Splits)"
                    )
                
                grid_seach_score = grid_result.best_score_
                print( "Optimized CV:", round( grid_seach_score, 3 ), "(MAE)")
                
                cv_output = pd.DataFrame(grid_result.cv_results_)
                cv_iteration = [
                    model,
                    n_estimators,
                    learning_rate,
                    max_features,
                    min_samples_split,
                    float(cv_output['mean_test_neg_mean_absolute_error'][grid_result.best_index_]),
                    float(cv_output['mean_train_neg_mean_absolute_error'][grid_result.best_index_]),
                    float(cv_output['mean_test_r2'][grid_result.best_index_]),
                    float(cv_output['mean_train_r2'][grid_result.best_index_])
                ]
                cv_archive.append( cv_iteration )
            
            # Identify "best-performing" machine-learning regressor (post-grid search)
            cv_summary = pd.DataFrame(cv_archive)
            cv_summary = cv_summary.rename(
                columns={
                    0: 'model',
                    1: 'n_estimators',
                    2: 'learning_rate',
                    3: 'max_features',
                    4: 'min_samples_split',
                    5: 'mae_test',
                    6: 'mae_train',
                    7: 'r2_test',
                    8: 'r2_train'
                }
            )
            # Sort models based on computed mean absolute error (MAE) metrics for the TRAIN and TEST matrices
            # Prioritize models with an MAE <= 0.5 for both the TRAIN/TEST matrices
            # However, if no model has an MAE <- 0.5 for the TRAIN matrix, then prioritize the model with the lowest MAE value (TRAIN)
            cv_summary['mae_train_abs'] = cv_summary['mae_train'].abs()
            cv_summary['mae_test_abs'] = cv_summary['mae_test'].abs()
            if cv_summary[ cv_summary['mae_train_abs'] <= 0.5 ].shape[0] == cv_summary.shape[0]:
                cv_summary = cv_summary.sort_values( by = ['mae_test_abs'], ascending = True ).reset_index( drop = True )
            elif cv_summary[ cv_summary[ 'mae_train_abs' ] <= 0.5 ].shape[0] > 0:
                cv_summary = cv_summary[ cv_summary[ 'mae_train_abs' ] <= 0.5 ]
                cv_summary = cv_summary.sort_values( by = ['mae_test_abs'], ascending = True ).reset_index( drop = True )
            else:
                cv_summary = cv_summary.sort_values( by = ['mae_train_abs'], ascending = True ).reset_index( drop = True )
            
            # Extract "optimal" model based on filtered/sorted 'cv_summary' dataframe
            opt_model = cv_summary.iloc[0]
            model = opt_model['model']
            
            # Apply optimized parameters/hyperparameters in expanded, bootstrapped model
            if str(model) == 'AdaBoostRegressor()':
                opt_model = AdaBoostRegressor(
                    n_estimators = int(opt_model['n_estimators']),
                    learning_rate = float(opt_model['learning_rate']),
                )
            elif str(model) == 'GradientBoostingRegressor()':
                opt_model = GradientBoostingRegressor(
                    n_estimators = int(opt_model['n_estimators']),
                    learning_rate = float(opt_model['learning_rate']),
                )
            elif str(model) == 'ExtraTreesRegressor()':
                opt_model = ExtraTreesRegressor(
                    n_estimators = int(opt_model['n_estimators']),
                    max_features = int(opt_model['max_features']),
                    min_samples_split = int(opt_model['min_samples_split']),
                    n_jobs = -1,
                )
            
            # Number of "fold" in cross-validation set by 'n_splits' parameter [Default: 5]
            # Number of times cross-validation is repeated set by 'n_repeats' parameter [Default: 10]
            print( "Initiate Cross Validation ...")
            cross_validation = RepeatedKFold(
                n_splits = 5,
                n_repeats = 200,
                random_state = 12345678
            )
            
            #cv_scores = cross_val_score(
            cv_scores = cross_validate(
                opt_model,
                X,
                Y,
                cv = cross_validation,
                n_jobs = 5,
                scoring = ['neg_mean_absolute_error','r2'],
                verbose = 1
            )
            
            print(
                'Model CV:',
                '[R2]',
                round( np.mean(cv_scores['test_r2']), 3 ),
                '(Mean) and',
                round( np.median(cv_scores['test_r2']), 3 ),
                '(Median)',
                '/ [MAE]',
                round( np.mean(cv_scores['test_neg_mean_absolute_error']), 3 ),
                '(Mean) and',
                round( np.median(cv_scores['test_neg_mean_absolute_error']), 3 ),
                '(Median)'
            )
            
            fitted_model = opt_model.fit(
                X,
                Y
            )
            
            # Prepare "predicted" X-dimensional dataframe to visualize predicted dimensional landscape
            Y_Predict = pd.DataFrame( opt_model.predict( X_Predict ))
            XY_Predict = pd.concat(
                [
                    X_Predict,
                    Y_Predict
                ],
                axis = 1
            )

            # Prepare "predicted" Z-coordinate dataframe for all line-items (e.g., HA variants by 'variant_hash')
            Z_Predict = pd.DataFrame( opt_model.predict( pXY ))
            XYZ_Predict = pd.concat(
                [
                    pXY,
                    Z_Predict
                ],
                axis = 1
            )
            
            # Reset column index: https://datascientyst.com/reset-column-names-index-pandas/
            XY_Predict.columns = range(XY_Predict.columns.size)
            XY_Predict = XY_Predict.rename(
                columns={
                    0: 'x',
                    1: 'y',
                    2: 'z'
                }
            )

            XYZ_Predict.columns = range(XYZ_Predict.columns.size)
            XYZ_Predict = XYZ_Predict.rename(
                columns={
                    0: 'x',
                    1: 'y',
                    2: 'z_predicted'
                }
            )
            
            # Prepare consolidated dataframe to export "predicted" 3D landscape
            # Establish length of index based on the number of predicted XYZ coordinates
            pred_landscape_iteration['x_predicted'] = XY_Predict['x']
            pred_landscape_iteration['y_predicted'] = XY_Predict['y']
            pred_landscape_iteration['z_predicted'] = XY_Predict['z']
            
            pred_landscape_iteration['subtype'] = str(composite.subtype.unique()[0])
            pred_landscape_iteration['protein'] = str(composite.protein.unique()[0])
            pred_landscape_iteration['x_factor'] = str(xfactors[i])
            pred_landscape_iteration['y_factor'] = str(yfactors[j])
            pred_landscape_iteration['z_factor'] = str(zfactors[k])
            pred_landscape_iteration['model'] = str(model)
            pred_landscape_iteration['n_estimators'] = int(cv_summary['n_estimators'][0])
            
            if pd.isnull( cv_summary['learning_rate'][0] ):
                pred_landscape_iteration['learning_rate'] = None
            else:
                pred_landscape_iteration['learning_rate'] = float(cv_summary['learning_rate'][0])
            if pd.isnull( cv_summary['max_features'][0] ):
                pred_landscape_iteration['max_features'] = None
            else:
                pred_landscape_iteration['max_features'] = int(cv_summary['max_features'][0])
            if pd.isnull( cv_summary['min_samples_split'][0] ):
                pred_landscape_iteration['min_samples_split'] = None
            else:
                pred_landscape_iteration['min_samples_split'] = int(cv_summary['min_samples_split'][0])
            
            pred_landscape_iteration['mae_test'] = float(cv_summary['mae_test'][0])
            pred_landscape_iteration['mae_train'] = float(cv_summary['mae_train'][0])
            pred_landscape_iteration['r2_test'] = float(cv_summary['r2_test'][0])
            pred_landscape_iteration['r2_train'] = float(cv_summary['r2_train'][0])
            pred_landscape_iteration['cv_mae_mean'] = float(round( np.mean(cv_scores['test_neg_mean_absolute_error']), 3 ))
            pred_landscape_iteration['cv_mae_median'] = float(round( np.median(cv_scores['test_neg_mean_absolute_error']), 3 ))
            pred_landscape_iteration['cv_r2_mean'] = float(round( np.mean(cv_scores['test_r2']), 3 ))
            pred_landscape_iteration['cv_r2_median'] = float(round( np.median(cv_scores['test_r2']), 3 ))

            # Prepare consolidated dataframe to export "predicted" Z-axis coordinates (reponse values) for applied XY model
            pred_metrics = pcomposite.copy()
            pred_metrics = pred_metrics[metadata[0:5]]
            pred_metrics['x_factor'] = str(xfactors[i])
            pred_metrics['y_factor'] = str(yfactors[j])
            pred_metrics['z_factor'] = str(zfactors[k])

            pred_metrics_output = pd.concat(
                [
                    pred_metrics,
                    XYZ_Predict
                ],
                axis = 1
            )

            # Export predicted XYZ landscapes to user-defined output TSV file
            if os.path.exists(outfile):
                with open(outfile,'a') as file:
                    pred_landscape_iteration.to_csv(
                    file,
                    sep='\t',
                    header=False,
                    index=False
                )
            else:
                with open(outfile,'w') as file:
                    pred_landscape_iteration.to_csv(
                    file,
                    sep='\t',
                    index=False
                )
                    
            # Export predicted Z-coordinate metric(s) to user-defined output TSV file
            if os.path.exists(pred_outfile):
                with open(pred_outfile,'a') as file:
                    pred_metrics_output.to_csv(
                    file,
                    sep='\t',
                    header=False,
                    index=False
                )
            else:
                with open(pred_outfile,'w') as file:
                    pred_metrics_output.to_csv(
                    file,
                    sep='\t',
                    index=False
                )
                
            # Summarize iteration ML analytics and output to TSV file
            pred_landscape_summary = pred_landscape_iteration.drop( ['x_predicted', 'y_predicted', 'z_predicted'], axis = 1 ).drop_duplicates().reset_index( drop = True )

            if os.path.exists(summary_outfile):
                with open(summary_outfile,'a') as file:
                    pred_landscape_summary.to_csv(
                    file,
                    sep='\t',
                    header=False,
                    index=False
                )
            else:
                with open(summary_outfile,'w') as file:
                    pred_landscape_summary.to_csv(
                    file,
                    sep='\t',
                    index=False
                )
