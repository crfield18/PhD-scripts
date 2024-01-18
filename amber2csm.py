'''
Re-format pdb files from MD simulations using AMBER and GLYCAM force fields via cpptraj.
Rename 3-letter sugar residue codes to match PDB for compatibility with CSM-carbohydrate,
a webserver predicting protein-carbohydrate binding affinity.
(Nguyen et al. 2022; https://biosig.lab.uq.edu.au/csm_carbohydrate/)

Args:
    --input/-i: Path to one or more pdb files

Returns:
  Formatted PDB file(s) are stored in a new or existing subdirectory named "amber2csm_output"
  located at the same path as the input PDB file.
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

    parser = ArgumentParser(description='Reformat pdb files for CSM-carbohydrate compatibility.')
    parser.add_argument('-i', '--input', type=Path, required=True, nargs='+',
                        help='Path to one or more pdb files')
    return parser.parse_args()

def line_filter(line:str):
    '''
    Format a line from a PDB file to be compatible with CSM-carbohydrate.
    beta-Glucose residues are given the residue code BGC and are moved to chain B.
    Protein residues are moved to chain A.

    Args:
        line:   String containing 1 line of an open PDB file.

    Returns:
        Formatted line string.
    '''

    amber_glc_codes = ['ROH', '4GA', '0GB']
    if line.startswith('ATOM'):
        if line[17:20] not in amber_glc_codes and line[21] == ' ':
            # Add Chain ID A to all non-glucose residues
            return f'{line[0:21]}A {line[23:-1]}\n'
        if line[17:20] in amber_glc_codes:
            # Replace GLYCAM glucose 3 letter codes with BGC and move to chain B
            return f'{line[0:17]}BGC B {line[23:-1]}\n'
        return line
    if line.startswith('TER'):
        return 'TER\n'
    return line

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
            output_file.writelines(line_filter(line))

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
