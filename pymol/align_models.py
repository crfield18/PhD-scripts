import glob
import json
import datetime
from argparse import ArgumentParser
from itertools import combinations
from math import isnan
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pymol import cmd
import psico.fullinit # Needed to add tmalign to PyMOL
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
# Look into radialtree to draw circular dendrograms


class PyMOLalign():
    def __init__(self, rms_results) -> None:
        # Raw RMSD results
        self.rms_results = rms_results
        # RMSD after refinement
        self.rmsd_refined = rms_results[0]
        # Number of aligned atoms after refinement
        self.nalign_atom_refined = rms_results[1]
        # Number of refinement cycles
        self.refine_cycles = rms_results[2]
        # RMSD before refinement
        self.rmsd_initial = rms_results[3]
        # Number of aligned atoms before refinement
        self.nalign_atom_initial = rms_results[4]
        # Raw alignment score
        self.score_initial = rms_results[5]
        # Number of residues aligned
        self.nalign_res = rms_results[6]

    def show_results(self):
        print(f'''
RMSD after refinement:\t\t\t{self.rmsd_refined}
# of aligned atoms after refinement:\t{self.nalign_atom_refined}
# of refinement cycles:\t\t\t{self.refine_cycles}
RMSD before refinement:\t\t\t{self.rmsd_initial}
# of aligned atoms before refinement:\t{self.nalign_atom_initial}
Raw alignment score:\t\t\t{self.score_initial}
# of residues aligned:\t\t\t{self.nalign_res}
''')

    def store_results(self, algorithm:str):
        return {algorithm:[self.rmsd_refined,self.nalign_res]}

class PyMOLcealign():
    def __init__(self, ce_results) -> None:
        # Raw RMSD results
        self.ce_results = ce_results
        # RMSD
        self.rmsd = self.ce_results['RMSD']
        # Alignment length
        self.nalign = self.ce_results['alignment_length']
        # Rotation matrix
        self.rotation_matrix = np.array(self.ce_results['rotation_matrix']).reshape(4,4)

    def show_results(self):
        print(f'''
RMSD after refinement:\t\t\t{self.rmsd}
# of aligned atoms after refinement:\t{self.nalign}
Model Rotation Matrix:\n{self.rotation_matrix}
''')

    def store_results(self):
        return {'cealign':[self.rmsd,self.nalign]}

class PyMOLAlignAll():
    def __init__(self, reference_object:str) -> None:
        self.ref = reference_object
        self.obj_list = self.non_sele_objs()
        self.results_dict = {}
        self.results_df = None

    def non_sele_objs(self):
        # Get a list of all object loaded into PyMOL, excluding 'sele'
        all_objects = cmd.get_object_list()
        non_sele_obj_list = [obj for obj in all_objects if obj != 'sele']
        return non_sele_obj_list

    def tmalign_all(self):
        for obj in self.obj_list:
            if obj != self.ref:
                tm = cmd.tmalign(obj, self.ref)
                sub_results_dict = {'tmalign': tm}
                self.results_dict.setdefault(obj, {}).update(sub_results_dict)

    def align_all(self):
        for obj in self.obj_list:
            if obj != self.ref:
                print(f'\nAligning based on sequence (align):\t{obj}')
                align_result_parser = PyMOLalign(cmd.align(obj, self.ref, cycles=10))
                align_result_parser.show_results()
                sub_results_dict = align_result_parser.store_results(algorithm='align')
                self.results_dict.setdefault(obj, {}).update(sub_results_dict)

    def cealign_all(self):
        for obj in self.obj_list:
            if obj != self.ref:
                print(f'\nAligning based on structure (cealign):\t{obj}')
                cealign_result_parser = PyMOLcealign(cmd.cealign(self.ref, obj))
                cealign_result_parser.show_results()
                sub_results_dict = cealign_result_parser.store_results()
                self.results_dict.setdefault(obj, {}).update(sub_results_dict)

    def super_all(self):
        for obj in self.obj_list:
            if obj != self.ref:
                print(f'\nSuperimposing (super):\t{obj}')
                align_result_parser = PyMOLalign(cmd.super(obj, self.ref))
                align_result_parser.show_results()
                sub_results_dict = align_result_parser.store_results(algorithm='super')
                self.results_dict.setdefault(obj, {}).update(sub_results_dict)

    def show_results_dict(self):
        print(self.results_dict)

    def get_results_dict(self):
        return self.results_dict

    def results_dict_to_df(self):
        # Convert the alignment results dictionary to a dataframe
        if not self.results_df:
            self.results_df = pd.DataFrame.from_dict(self.results_dict, orient='index')
            # Expand lists and assign elements to new columns
            self.results_df[['align_RMSD',
                             'align_nres']] = self.results_df['align'].apply(pd.Series)
            self.results_df[['cealign_RMSD',
                             'cealign_nres']] = self.results_df['cealign'].apply(pd.Series)
            self.results_df[['super_RMSD',
                             'super_nres']] = self.results_df['super'].apply(pd.Series)

            # Drop the original columns with lists
            self.results_df.rename(columns={'tmalign': 'TM-align'}, inplace=True)
            self.results_df.drop(['align', 'cealign', 'super'], axis=1, inplace=True)

        return self.results_df

class TMalignMatrix():
    def __init__(self, align_all_results:dict) -> None:
        self.results = align_all_results
        self.tmmatrix = pd.DataFrame()

    def get_sequence_length(self, object_name):
        # Count the number of CA atoms in given object
        sequence_length = cmd.count_atoms(f'{object_name} and name CA')
        return sequence_length

    def calculate_matrix(self, save_temp:bool):
        # Generate a matrix of TMalign scores to build a structural dendrogram
        # Create a list of unique PDB codes ordered by their TMalign score
        model_list = [model for model, results in self.results.items()]
        unique_model_list = get_unique_models(model_list, top_x=100)

        try:
            self.tmmatrix = pd.read_csv(glob.glob("*tmalign-matrix*.csv")[0], index_col=0)
        except IndexError:
            self.tmmatrix = pd.DataFrame(index=unique_model_list, columns=unique_model_list)

        self.show_matrix()

        # Calculate TMalign scores for all unique combinations of models
            # This would be more efficient if storing values in a list/dict
            # Directly updating the dataframe makes the data easier to parse visually
            # The real bottleneck here is going to TMalign anyway so I'm leaving it
        for pair in combinations(unique_model_list, 2):
            print(pair)
            if isnan(self.tmmatrix.loc[pair[0], pair[1]]) or isnan(self.tmmatrix.loc[pair[1], pair[0]]):
                # Set the smaller object as the reference (target) object
                object_1_length = self.get_sequence_length(pair[0])
                object_2_length = self.get_sequence_length(pair[1])

                if object_1_length >= object_2_length:
                    target = pair[1]
                    mobile = pair[0]
                elif object_1_length < object_2_length:
                    target = pair[0]
                    mobile = pair[1]
                
                tm = cmd.tmalign(mobile, target)

                self.tmmatrix.loc[pair[0], pair[1]] = tm
                self.tmmatrix.loc[pair[1], pair[0]] = tm

                self.show_matrix()
    
                # Save results as they are calculated to cut down on recalculating TMalign scores
                # after a crash or a restart
                if save_temp:
                    df_to_file(self.tmmatrix, filename='tmalign-matrix', filetype='csv')

        # Fill the diagonal with TMscores of 1 (always the score for identical proteins)
            # ∴ We do not need to run any calculations
        for obj in unique_model_list:
            self.tmmatrix.loc[obj,obj] = 1.0
        
        self.show_matrix()
        return self.tmmatrix
    
    def show_matrix(self):
        print(self.tmmatrix)

    def get_matrix(self):
        return self.tmmatrix
    
    def make_dendrogram(self):
        # Invert the DataFrame (subtract each value from 1)
        inverted_df = 1 - self.tmmatrix

        # Convert full distance matrix to condensed matrix
        similarity_matrix = inverted_df.values
        condensed_matrix = squareform(similarity_matrix)

        # Compute hierarchical clustering
        # average = UPGMA
        # single = nearest neighbour/minimum evolution (good for large datasets)
        Z = linkage(condensed_matrix, method='average')

        # Plot dendrogram
        plt.figure(figsize=(8, 6))
        dendrogram(Z, orientation='left', labels=inverted_df.index.tolist(), leaf_font_size=8)
        plt.title('TMalign Score (Structural) Dendrogram')
        axes = plt.gca()
        axes.set_xlim(1,0)
        # plt.show()
        plt.savefig('tmalign-score-dendrogram.png', transparent=None, dpi=300)

def get_unique_models(model_list:list, top_x:int=25):
    # Used for checking duplicates
    unique_model_set = set()
    unique_model_list = []

    # Create a list of unique PDB codes ordered by their TMalign score
    for m in model_list:
        if len(unique_model_list) < top_x:
            # Could check for duplicates of first 2 letters of pdb code to find proteins from different sources
            if m.split('-')[0] not in unique_model_set:
                unique_model_set.add(m.split('-')[0])
                unique_model_list.append(m)
        else:
            break
 
    return unique_model_list

def df_to_file(dataframe:pd.DataFrame, filename:str, filetype:str='csv', index_label:str=None):
    # Write alignment results dataframe to a csv (default) or excel file
    method_map = {
        'csv': dataframe.to_csv,
        'xlsx': dataframe.to_excel
    }
    
    current_date = datetime.datetime.now().strftime('%Y%m%d')
    if filetype not in ['csv', 'xlsx']:
        raise ValueError('Invalid filetype specified. Please use "csv" or "xlsx".')
    if index_label:
        method_map[filetype](f'{filename}-{current_date}.{filetype}', index=True, index_label=index_label)
    else:
        method_map[filetype](f'{filename}-{current_date}.{filetype}', index=True)

def script_args():
    parser = ArgumentParser()
    parser.add_argument('-r', '--reference', type=Path, required=True,
                        help='Path to reference pdb/cif file.')
    parser.add_argument('-d', '--directory', type=Path, required=True,
                        help='Path to directory containing all pdb/cif files to align to the reference structure.')
    return parser.parse_args()

def align():
    user_args = script_args()

    # Get a list of all the files in the user specified directory
    test_dir_path = user_args.directory.resolve()
    cif_list = [file for file in test_dir_path.iterdir()
                if file.is_file() and file.suffix in ('.pdb', '.cif')]

    # Get the path for the user specified reference pdb/cif file
    cif_ref_path = user_args.reference.resolve()

    # Load all pdb/cif files in cwd
    if cif_ref_path not in cif_list:
        cmd.load(cif_ref_path)

    for cif in cif_list:
        print(f'Loading: {cif}')
        cmd.load(cif)

    # Remove waters
    cmd.remove('resn hoh')

    # Initialise a PyMOL instance and set the reference object
    pymol_instance = PyMOLAlignAll(cif_ref_path.stem)

    # Align all other objects to the reference using align, cealign, super and tmalign
    for alignment_algorithm in (pymol_instance.align_all,
                                pymol_instance.cealign_all,
                                pymol_instance.super_all,
                                pymol_instance.tmalign_all):
        try:
            print('-' * 100)
            alignment_algorithm()
        except:
            pass

    # Sort results_dict by TMalign score
    pymol_instance.results_dict = dict(sorted(pymol_instance.results_dict.items(), key=lambda item: item[1]['tmalign'], reverse=True))

    # Write results to a JSON file (mainly for testing)
    with open('results.json', 'w', encoding='UTF8') as json_file:
        json.dump(pymol_instance.results_dict, json_file, indent=4)

    # Convert the results dictionary to a dataframe and write to csv and xlsx files
    pymol_instance.results_dict_to_df()
    for extension in ('csv', 'xlsx'):
        df_to_file(dataframe=pymol_instance.results_df,
                                filename=f'{pymol_instance.ref}-alignment',
                                index_label='PDB_code_with_chain',
                                filetype=extension)


    cmd.reinitialize()

    # Calculate TMalign score matrix
        # Used to create a structural dendrogram

    tmmatrix = TMalignMatrix(pymol_instance.results_dict)

    tmmatrix.calculate_matrix(save_temp=True)

    for extension in ('csv', 'xlsx'):
        df_to_file(dataframe=tmmatrix.get_matrix(),
                                filename='tmalign-matrix',
                                filetype=extension)

    tmmatrix.make_dendrogram()

    return pymol_instance.ref, cif_ref_path, pymol_instance.results_dict

def visualise_top_x(results:dict, top_x:int, reference_obj:str):
    # Create PyMOL sessions for the top x PDBs with the highest TMalign scores
    # Only 1 chain per PDB code is used (often there are identical chains in the same PDB) 

    # Create a list of unique PDB codes ordered by their TMalign score
    model_list = [model for model, results in results.items()]
    unique_model_list = get_unique_models(model_list, top_x)

    # Open and align copies of the chain for each alignment method
    # PLACEHOLDER! DOES NOT WORK YET!
    for m in unique_model_list:
        cmd.load(reference_obj)
        # CHANGE: Name objects '{algorithm}_{rmsd/tm to 4 d.p.}
        for method in ('align', 'cealign', 'super', 'tmalign'):
            cmd.load(f'{m}.pdb', f'{m}_{method}_{m}')

        # print(cmd.get_names())
        cmd.remove('resn hoh')
        cmd.cealign(reference_obj.stem, f'{m}_cealign')
        cmd.align(f'{m}_align', reference_obj.stem)
        cmd.super(f'{m}_super', reference_obj.stem)
        cmd.tmalign(f'{m}_tmalign', reference_obj.stem)

        # 0.65 = 35% visible
        cmd.set('cartoon_transparency', 0.65, f'(all and not {reference_obj.stem})')
        cmd.center(reference_obj.stem)
        
        # Save pymol session file
        cmd.save(f'{m}.pse')
        cmd.reinitialize()

def main():
    reference_obj, reference_obj_path, alignment_results = align()

    visualise_top_x(alignment_results, 25, reference_obj_path)

if __name__ == '__main__':
    main()

# Jobs for Saturday:
    # Separate results generation to its own function probably
    # Add reminder sheet to excel doc with basic info on each alignment algorithm/how to interpret scores
    # name objects in visualise_top_x based on RMSD/TMalign score
