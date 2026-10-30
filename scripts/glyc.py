#!/usr/bin/env python3
import warnings
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, Selection, DSSP
from Bio.PDB.DSSP import residue_max_acc
import pandas as pd
from sklearn.neighbors import radius_neighbors_graph
from pathlib import Path
import argparse


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_pdb_file", required=True,help="Input PDB file")

    parser.add_argument("-v", "--variant_hash", required=True,
        help="Alphanumeric character STRING linked to a specific amino acid sequence")

    parser.add_argument("-t", "--timestamp_version", required=True,
        help="System TIMESTAMP applied to distinct analytical applications")

    parser.add_argument("-o", "--output_csv_file", required=True,
        help="""Output csv file. If no argument given, then csv file with the same name as 
              the input is written to same directory""")

    return parser


def ignore_pdb_warnings():
    """
    Supress the biopython pdb import warnings
    """
    warnings.simplefilter("ignore", BiopythonWarning)


def add_node_id_column(df):
    """
    Add a node id column. A node id is a unique identifier for each residue comprised of
    chain, residue number, intertion code, and residue name
    """
    return df["chain"] + ":" + df["site"].astype(str) + ":" + \
                                df["insertion_code"] + ":" + df["residue"]


def make_PDB_df(structure):
    """
    Create a dataframe of residue properties contained with the PDB strucure input
    """
    structure_list = []
    for atom in Selection.unfold_entities(structure, "A"):
        atom_name = str(atom.id)
        if atom_name == "CA":
            resnum = str(atom.get_parent().id[1])
            resname = str(atom.get_parent().resname)
            bfactor = float(atom.bfactor)
            if atom.get_parent().id[2] ==  ":":
                insertion_code = " "
            else:
                insertion_code = atom.get_parent().id[2]
            chain = str(atom.get_parent().get_parent().id)
            node_id = chain + ":" + resnum + ":" + insertion_code + ":" + resname
            structure_list.append([node_id, atom.coord[0], atom.coord[1], atom.coord[2], bfactor])

    structure_df = pd.DataFrame(structure_list, columns=["node_id", "x_coordinate", "y_coordinate", "z_coordinate", "b_factor"])
    structure_df["chain"] = structure_df.node_id.str.split(":", expand=True)[0]
    structure_df["site"] = structure_df.node_id.str.split(":", expand=True)[1].astype(int)
    structure_df["insertion_code"] = structure_df.node_id.str.split(":", expand=True)[2]
    structure_df["residue"] = structure_df.node_id.str.split(":", expand=True)[3]

    return structure_df


def add_glycosylation_motif_col(df, column):
    """
    Add a column called "glycosylation_motif" with values NXS or NXP as to whether the residue is 
    in the ASN-X(not Pro)-S/T of the putative glycosylation motif. Null if otherwise
    """

    df.loc[(df[column] == "ASN") & 
           (df[column].shift(-1) != "PRO") & 
           (df[column].shift(-2) == "SER"), "glycosylation_motif"] = "NXS"
    df.loc[(df[column] != "PRO") & 
           (df[column].shift(1) == "ASN") & 
           (df[column].shift(-1) == "SER"), "glycosylation_motif"] = "NXS"
    df.loc[(df[column] == "SER") & 
           (df[column].shift(1) != "PRO") & 
           (df[column].shift(2) == "ASN"), "glycosylation_motif"] = "NXS"

    df.loc[(df[column] == "ASN") & 
           (df[column].shift(-1) != "PRO") & 
           (df[column].shift(-2) == "THR"), "glycosylation_motif"] = "NXT"
    df.loc[(df[column] != "PRO") & 
           (df[column].shift(1) == "ASN") & 
           (df[column].shift(-1) == "THR"), "glycosylation_motif"] = "NXT"
    df.loc[(df[column] == "THR") & 
           (df[column].shift(1) != "PRO") & 
           (df[column].shift(2) == "ASN"), "glycosylation_motif"] = "NXT"



def calculate_distance_to_glycosylation_motif(df, glyc_motif_col='glycosylation_motif', radius=20):
    """
    For every glycosylation motif in the structure/dataframe, add a column that is the distance from
    the residue of the column to that particular Asn of the putative glycosylation motif. Columns will 
    look like "Distance_to_N<residue_number>:<chain id>". Also adds a column "glycosylation site" which 
    is the residue number of the Asn of the putative glycosylation site that is closest to the residue
    of the row, and the "glycosylation_distance" column is how far away aforementioned Asn is from the residue
    """
    model_coords = df[['x_coordinate', 'y_coordinate', 'z_coordinate']].values
    glyc_sites = df.loc[(~df[glyc_motif_col].isna()) & (df.residue == 'ASN')].copy()
    glyc_sites['GlySite'] = 'Distance_to_N' + glyc_sites['site'].astype(str) + ':' + glyc_sites['chain']
    cols_rename_dict = glyc_sites['GlySite'].to_dict()
    m = radius_neighbors_graph(model_coords, radius, mode='distance', include_self=True)
    distance_matrix = pd.DataFrame(m.toarray()).replace(0, radius)
    distance_matrix.rename(columns=cols_rename_dict, inplace=True)
    distance_matrix = distance_matrix[list(cols_rename_dict.values())].copy()
    df = df.merge(distance_matrix, left_index=True, right_index=True)
    glycosylation_cols = [col for col in df if 'Distance_to_N' in col]
    for glyc_motif_col_name in glycosylation_cols:
        df.loc[(df.site == int(glyc_motif_col_name.split(':')[0][13:])) & (df.chain == glyc_motif_col_name[-1]), glyc_motif_col_name] = 0.0
    df['glycosylation_site'] = df[glycosylation_cols].idxmin(axis=1)
    df['glycosylation_site'] = df['glycosylation_site'].str.replace(r'Distance_to_N(\d+):\D', r'\1', regex=True)
    df['glycosylation_distance'] = df[glycosylation_cols].min(axis=1)
    df.loc[df['glycosylation_distance'] == 20.0, 'glycosylation_site'] = None
    df.loc[df['glycosylation_distance'] == 20.0, 'glycosylation_distance'] = None

    return df


def add_surface_accessibility_column(df):
    """
    Add a column for whether or not the residue is exposed or buried based on whether its
    rsa is above/below 0.25 and is more than 20 Angstroms from a glycosylation site
    """
    df['surface_accessibility'] = ['accessible' if x >= 0.25 else 'inaccessible' for x in df["rsa"]]
    df['glycan_shielding'] = ['exposed' if x >= 19 else 'glycan_shielded' for x in df["glycosylation_distance"]]
    
    
def add_disulfide_bond_column(df):
    """
    Placeholder column, probably will have this information in another table
    """
    df["disulfide_bond"] = False


def remove_unnecessary_glyc_columns(df):
    """
    Removes columns related to glycosylation distance that are needed in the final output. They could be useful for
    future purposes though.
    """
    unnecessary_glyc_columns = [col for col in df.columns if "Distance_to" in col]
    df.drop(columns=unnecessary_glyc_columns, inplace=True)


def create_dssp_df(structure, file, columns_to_keep=['secondary_structure', 'asa', 'rsa']):
    """
    Creates a dataframe from the output of the program DSSP. These calculations are secondary structure, 
    solvent accessibility, and relative solvent accessibility for each residue. However, more data exists 
    in the DSSP output.

    structure           = PDB structure parsed with biopython
    file                = the pdb file that "structure" was parsed from
    columns_to_keep     = columns from DSSP output to keep. Available columns:
        chain, HetFlag, site, insertion_code, Index, residue, 
        secondary_structure (renamed from SS), rsa, Phi, Psi, 
        NH-->O_1_relidx, NH-->O_1_energy, O-->NH_1_relidx, O-->NH_1_energy, 
        NH-->O_2_relidx, NH-->O_2_energy, O-->NH_2_relidx, O-->NH_2_energy
        asa is also an option to keep, this is a column thats calculated after DSSP

    """
    dssp_data = DSSP(structure[0], file, dssp='mkdssp')
    dssp_df = pd.DataFrame([tuple(x[0]) + tuple(x[1]) + y for x, y in list(zip(dssp_data.keys(), list(dssp_data)))], 
        columns=['chain', 'HetFlag', 'site', 'insertion_code','Index', 'residue', 'SS', 'rsa', 'Phi', 'Psi', 
                'NH-->O_1_relidx', 'NH-->O_1_energy','O-->NH_1_relidx', 'O-->NH_1_energy', 
                'NH-->O_2_relidx', 'NH-->O_2_energy','O-->NH_2_relidx', 'O-->NH_2_energy'])
    aa_mapper = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS', 
                'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
                'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
                'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}
    insertioncode_mapper = {' ': ' '}
    dssp_df = dssp_df.replace({"insertion_code": insertioncode_mapper})
    dssp_df = dssp_df.replace({"residue": aa_mapper})
    dssp_df.rename(columns={'SS': 'secondary_structure'}, inplace=True)
    dssp_df['node_id'] = add_node_id_column(dssp_df)
    dssp_df['maxasa'] = dssp_df['residue'].map(residue_max_acc['Sander'])
    dssp_df['asa'] = (dssp_df['maxasa'] * dssp_df['rsa']).astype(int)
    dssp_df['site'] = dssp_df['site'].astype('int')
    dssp_df = dssp_df[['node_id'] + columns_to_keep]
    #dssp_df = dssp_df['secondary_structure','asa','rsa','maxasa','site','residue','insertion_code','secondary_structure','node_id']

    return dssp_df


def format_df_for_output(df, column_order=["variant_hash","chain","site","residue","insertion_code","node_id",
                                            "x_coordinate","y_coordinate","z_coordinate",
                                            "secondary_structure","disulfide_bond",
                                            "glycosylation_site","glycosylation_motif",
                                            "surface_accessibility","glycan_shielding","asa","rsa","b_factor",
                                            "glycosylation_distance","version"]):
    """
    Format the df to match the SQL schema
    """
    #float32_cols = list(df.select_dtypes(include='float32'))
    #df[float32_cols] = df[float32_cols].astype('float64')


    df = df[column_order]

    return df


def main():

    args = get_parser().parse_args()
    pdb_file            = Path(args.input_pdb_file)
    variant_hash        = args.variant_hash
    timestamp_version   = args.timestamp_version
    if args.output_csv_file:
        output_csv_file = Path(args.output_csv_file)
    else:
        output_csv_file = Path(pdb_file.parent.joinpath(str(pdb_file)+"glyc_ASA_dist.tsv"))

    ignore_pdb_warnings()
    parser = PDBParser()
    structure = parser.get_structure(pdb_file.stem, pdb_file)
    structure_df = make_PDB_df(structure)
    print(structure_df)
    add_glycosylation_motif_col(structure_df, "residue")
    print("add residue")
    print(structure_df)
    add_disulfide_bond_column(structure_df)
    print("add disulfide")
    print(structure_df)
    structure_df = calculate_distance_to_glycosylation_motif(structure_df)
    print("add glyc dist")
    print(structure_df)
    remove_unnecessary_glyc_columns(structure_df)
    print("remove Distance tos")
    print(structure_df)
    dssp_df = create_dssp_df(structure, pdb_file)
    print("dssp dfs")
    print(dssp_df)
    #structure_df = pd.merge(structure_df,dssp_df, on="node_id", how="left",verify_integrity=True)
    structure_df = structure_df.merge(dssp_df,on="node_id",how="inner",sort=True,validate="1:1" )
    print("merge dfs")
    print(structure_df)
    add_surface_accessibility_column(structure_df)
    print("SA")
    print(structure_df)

    if variant_hash:
        structure_df["variant_hash"] = variant_hash
    if timestamp_version:
        structure_df["version"] = pd.to_datetime(timestamp_version)

    structure_df = format_df_for_output(structure_df)
    structure_df.to_csv(output_csv_file, sep ="\t",header=False,index=False)
    print(structure_df.dtypes)


if __name__ == "__main__":
    main()