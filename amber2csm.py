'''
Re-format pdb files from MD simulations using AMBER and GLYCAM force fields via cpptraj.
Rename 3-letter sugar residue codes to match PDB for compatibility with CSM-carbohydrate,
a webserver predicting protein-carbohydrate binding affinity.
(Nguyen et al. 2022; https://biosig.lab.uq.edu.au/csm_carbohydrate/)

Args:
    --input/-i: Path to one or more pdb files

Returns:
    "config.txt" in the input file's directory
'''

import os
from pathlib import Path
from argparse import ArgumentParser

def script_args():
    '''
    Parse command-line arguments and return a Namespace object containing
    pathlib.Path paths to each user-defined input file.

    Returns:
        argparse.Namespace: A Namespace object containing input file paths.
    '''
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', type=Path, required=True, nargs='+', help='')
    return parser.parse_args()

def convert_to_csm(pdb_file):
    '''
    Formats each user-defined PDB file to be compatible with CSM-carbohydrate.
    beta-Glucose residues are given the residue code BGC and are moved to chain B.

    Args:
        input_file:   pathlib.Path path to one of the input PDB files.

    Returns:
        None
    '''

    # Get directory with all the pdb files in
    pdb_parent_dir = pdb_file.resolve().parent

    # Make an output directory inside the pdb file directory
    results_dir = pdb_parent_dir / f'{os.path.basename(__file__).split(".")[0]}_output'
    results_dir.mkdir(parents=True, exist_ok=True)

    output_file_path = results_dir / f'{pdb_file.stem}_output.pdb'
    with open(pdb_file, 'r', encoding='utf-8') as input_file, \
         open(output_file_path, 'w', encoding='utf-8') as output_file:
        for line in input_file:
            amber_glc_codes = ['ROH', '4GA', '0GB']
            if line[0:3] == 'TER':
                output_file.writelines('TER\n')
            elif line[0:4] != 'ATOM' or line[17:20] not in amber_glc_codes:
                output_file.writelines(line)
            else:
                newline = f'{line[0:17]}BGC B {line[23:-1]}\n'
                # print(newline)
                output_file.writelines(newline)

def main():
    '''
    Primary execution point of the script.

    Usage:
        To run this script, execute it from the command line as follows:
        $ python3 amber2csm.py -i [path/to/file]
    '''

    user_args = script_args()
    for file in user_args.input:
        convert_to_csm(file)

if __name__ == '__main__':
    main()
