from pathlib import Path
from argparse import ArgumentParser
import json
# Figure out how to only use pandas
import numpy as np
import pandas as pd
from pymol import cmd
import psico.fullinit
import datetime
from math import isnan
from itertools import permutations

results_dict = {}

class PyMOLRMS():
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

    @property
    def rmsd(self):
        return self.rmsd_refined

class PyMOLcealign():
    def __init__(self, ce_results) -> None:
        # Raw RMSD results
        self.ce_results = ce_results
        # RMSD
        self.rmsd = self.ce_results['RMSD']
        # Alignment length
        self.nalign = self.ce_results['alignment_length']
        # Rotation matrix
        self.rotation_matrix = np.array(self.ce_results['rotation_matrix'])

class PyMOLAlign():
    def __init__(self, reference_object:str) -> None:
        self.ref = reference_object
        self.obj_list = PyMOLAlign.non_sele_objs()
        self.results_dict = {}

    def non_sele_objs():
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
                rms = cmd.align(obj, self.ref)
                align_result_parser = PyMOLRMS(rms)
                align_result_parser.show_results()
                sub_results_dict = {'align': [align_result_parser.rmsd_refined,
                                              align_result_parser.nalign_res]}
                self.results_dict.setdefault(obj, {}).update(sub_results_dict)

    def cealign_all(self):
        for obj in self.obj_list:
            if obj != self.ref:
                print(f'\nAligning based on structure (cealign):\t{obj}')
                rms = cmd.cealign(self.ref, obj)
                print(rms)
                sub_results_dict = {'cealign': [rms['RMSD'], rms['alignment_length']]}
                self.results_dict.setdefault(obj, {}).update(sub_results_dict)

    def super_all(self):
        for obj in self.obj_list:
            if obj != self.ref:
                print(f'\nSuperimposing (super):\t{obj}')
                align_result_parser = PyMOLRMS(rms_results=cmd.super(obj, self.ref))
                align_result_parser.show_results()
                sub_results_dict = {'super': [align_result_parser.rmsd_refined,
                                              align_result_parser.nalign_res]}
                self.results_dict.setdefault(obj, {}).update(sub_results_dict)

    def tm_matrix(self):


        # Used for checking duplicates
        unique_model_set = set()
        unique_model_list = []

        # Create a list of unique PDB codes ordered by their TMalign score
        model_list = [model for model, results in self.results_dict.items()]
        for m in model_list:
            if len(unique_model_list) < 25:
                # Could check for duplicates of first 2 letters of pdb code to find proteins from different sources
                if m.split('-')[0] not in unique_model_set:
                    unique_model_set.add(m.split('-')[0])
                    unique_model_list.append(m)
                else:
                    pass
            else:
                break

        print(unique_model_list)

        self.tmmatrix = pd.DataFrame(index=unique_model_list, columns=unique_model_list)
        for pair in permutations(unique_model_list, 2):
            print(pair)
            tm = cmd.tmalign(pair[0], pair[1])

            if tm > self.tmmatrix.loc[pair[0], pair[1]] or isnan(self.tmmatrix.loc[pair[0], pair[1]]):
                self.tmmatrix.loc[pair[0], pair[1]] = tm
                self.tmmatrix.loc[pair[1], pair[0]] = tm
            PyMOLAlign.tm_matrix_show(self)

        for obj in unique_model_list:
            self.tmmatrix.loc[obj,obj] = 1.0
        return self.tmmatrix
    
    def tm_matrix_show(self):
        print(self.tmmatrix)
    
    def tm_matrix_to_csv(self):
        # csv_file = 'tmalign_matrix.csv'
        self.tmmatrix.to_csv('tmalign_matrix.csv', index=True)
        # return csv_file
    

def script_args():
    parser = ArgumentParser()
    parser.add_argument('-r', '--reference', type=Path, required=True,
                        help='Path to reference pdb/cif file.')
    parser.add_argument('-d', '--directory', type=Path, required=True,
                        help='Path to directory containing all pdb/cif files to align to the reference structure.')
    return parser.parse_args()

def visualise_top_x(results:dict, top_x:int, reference_obj:str):
    # Used for checking duplicates
    unique_model_set = set()
    unique_model_list = []

    # Create a list of unique PDB codes ordered by their TMalign score
    model_list = [model for model, results in results.items()]
    for m in model_list:
        if len(unique_model_list) < top_x:
            # Could check for duplicates of first 2 letters of pdb code to find proteins from different sources
            if m.split('-')[0] not in unique_model_set:
                unique_model_set.add(m.split('-')[0])
                unique_model_list.append(m)
            else:
                pass
        else:
            break

    print(unique_model_list)

    # Open and align copies of the chain for each alignment method
    # PLACEHOLDER! DOES NOT WORK YET!
    for m in unique_model_list:
        cmd.load(f'{reference_obj}.pdb')
        # CHANGE: Name objects '{algorithm}_{rmsd/tm to 4 d.p.}
        for method in ('align', 'cealign', 'super', 'tmalign'):
            cmd.load(f'{m}.cif', f'{m}_{method}_{m}')

        # print(cmd.get_names())
        cmd.remove('resn hoh')
        cmd.cealign(reference_obj, f'{m}_cealign')
        cmd.align(f'{m}_align', reference_obj)
        cmd.super(f'{m}_super', reference_obj)
        cmd.tmalign(f'{m}_tmalign', reference_obj)


        # Set cartoon transparency to 0.65 (35% visible)
        cmd.set('cartoon_transparency', 0.65, f'(all and not {reference_obj})')

        cmd.center(reference_obj)
        # Save pymol session file
        cmd.save(f'{m}.pse')
        cmd.reinitialize()

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

    pymol_instance = PyMOLAlign(cif_ref_path.stem)

    for alignment_algorithm in (pymol_instance.align_all, pymol_instance.cealign_all, pymol_instance.super_all, pymol_instance.tmalign_all):
        print('-' * 100)
        alignment_algorithm()

    # Write results to a JSON file (mainly for testing)
    with open('results.json', 'w', encoding='UTF8') as json_file:
        results_sorted_by_tm = dict(sorted(pymol_instance.results_dict.items(),
                                           key=lambda item: item[1]["tmalign"], reverse=True))
        json.dump(results_sorted_by_tm, json_file, indent=4)

    pymol_instance.tm_matrix()
    pymol_instance.tm_matrix_to_csv()

    return pymol_instance.ref, pymol_instance.results_dict

if __name__ == '__main__':
    reference_obj, alignment_results = align()
    with open('results.json', 'r') as results_json:
        alignment_results = json.load(results_json)
    reference_obj = 'ranked_0'
    visualise_top_x(alignment_results, 25, reference_obj)

    # Write results to a csv file
    df = pd.DataFrame.from_dict(alignment_results, orient='index')

    print(df)
    # Expand lists and assign elements to new columns
    df[['align_RMSD', 'align_nres']] = df['align'].apply(pd.Series)
    df[['cealign_RMSD', 'cealign_nres']] = df['cealign'].apply(pd.Series)
    df[['super_RMSD', 'super_nres']] = df['super'].apply(pd.Series)

    # Drop the original columns with lists
    df.rename(columns={'tmalign': 'TM-align'}, inplace=True)
    df.drop(['align', 'cealign', 'super'], axis=1, inplace=True)

    # Print the DataFrame
    print(df)

    current_date = datetime.datetime.now().strftime('%Y%m%d')
    df.to_csv(f'{reference_obj}-alignment-{current_date}.csv', index=True, index_label='PDB_code_with_chain')
    df.to_excel(f'{reference_obj}-alignment-{current_date}.xlsx', index=True, index_label='PDB_code_with_chain')


# Jobs for saturday:
    # Lint Code
    # Separate results generation to its own function probably
    # Add reminder sheet to excel doc with basic info on each alignment algorithm/how to interpret scores