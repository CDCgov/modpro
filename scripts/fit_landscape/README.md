# DSU Fitness Landscapes: Visualizing Influenza Virus Evolutionary Trajectories in Multidimensional Landscapes

## Project Overview

The NCIRD Influenza Division (ID) is recognized globally for its leadership in influenza virus surveillance and prevention policy. Initiated by the “Sequence First” Initiative in 2016, the effort has greatly expanded the use of Next Generation Sequencing (NGS) to sequence all isolated specimens within the national surveillance network, significantly enhancing viral characterization efficiency. From 2018 to 2022, collaborative initiatives between CDC labs and National Influenza Reference Centers (NIRC) have processed between 4600 and 10900+ specimens annually, detecting 1000 – 3200+ unique HA protein viral variants, all validated across 6200 – 57000+ test results.

Despite comprehensive strategies to manage and integrate sequence surveillance data within the NCIRD Data Warehouse (NDW), transforming these complex datasets into accessible analytics and visuals to aid decision-making remains a significant challenge. Additionally, any developed resources need to be both dynamic and scalable to adapt to the rapidly changing influenza virus landscape.

### Problem Statement

The primary challenge is to simplify and visualize complex sequence surveillance data to assist in leadership decision-making and inform public health policies effectively.

### Solution Approach

The DSU Fitness Landscapes project employs fitness landscape modeling to depict complex virological data relationships. A 3D modeled surface, created by overlaying a two-dimensional plane of continuous data distributions with an independent "fitness" variable, is sliced into 2D plots along a summary path using locally-weighted regression. This method produces visualizations that map the fitness peaks and troughs of the data, thereby illustrating evolutionary "fitness" of various viral strains.

## Project Structure

- `bin/`: Contains executable scripts and binaries.
- `data/alphafold2/`: Houses raw and processed data files.
- `sql/`: SQL scripts for database management and queries.
- `tbl/`: Stores tables and other structured data files.
- `FitnessLandscape`: Bash wrapper
## Team Composition

- **Administrators**: Responsible for project oversight, management tasks, and coordination.
  - Brian Mann
  - Brian Lee
- **Developers**: Focused on developing, testing, and maintaining the software components of the project.
  - Stacy Carter
  - Akilia Mathie
  - Samuel Wiley
## Project Objectives

1. Model genetic and phenotypic evolution of the HA protein.
2. Predict evolutionary trajectories and potential pathways.
3. Direct specimen prioritization for antigenic surveillance.
4. Inform selection of candidate vaccine viruses (CVVs) to enhance responses against dominant influenza strains.
5. Utilize historical data from 2018-2022 and current data from the 2023/24 Northern Hemisphere influenza season to build and refine predictive models.

## Usage Information

The `Fitness Landscape` Linux shell BASH script includes several commands:

````{bash}
    Fitness Landscape   
	            autoML [INPUT] [OUTPUT] [DELIMITER] [COMPOSITE|MIXED] [ALL| 
                ANTIGENIC|PROTEOMIC]
                distance
                do
                impute
                merge
                predictML [INPUT] [OUTPUT] [DELIMITER] [ITERATIVE] 
                [ALL|ANTIGENIC|MIXED|PROTEOMIC] [REFERENCE INPUT]
                update               
                        
                        
````

| Subcommand | Description |
| :- | :- |
| `autoML` | Identifies optimal ML models for user-defined datasets which allows for regressor selection.|
| `distance` | Computes structural distance metric.|
| `do` | Performs sequential execution of SQL queries in the CDP-CLIENT-02 environment via the [**hadoopBySSH**](https://git.biotech.cdc.gov/flu-informatics/hadoopBySSH) `himpala` and BASH scripts.|
| `impute` | Imputes missing data in local files. |
| `merge` | Merges local data files. |
| `predictML` | Model predicts/generates optimized XYZ -coordinate landscapes using the regressors identified in the autoML.|
| `update` | Refreshes local data files. |
|  |  |