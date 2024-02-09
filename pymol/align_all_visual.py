'''
Align multiple PDB files using PyMOL

It takes input PDB files, aligns them to a reference structure (if provided),
and saves the PyMOL session as a .pse file.

Usage:
    python align_pdb.py -i input1.pdb input2.pdb -r reference.pdb

Args:
    -i, --input:    Input PDB files to align (≥2).
    -g, --gui:      Run PyMOL instances with the GUI.
    -r, --ref:      Reference PDB file to align to (1).
'''

from pathlib import Path
from argparse import ArgumentParser
import pymolPy3 as pymol

def script_args():
    '''
    Parse command-line arguments.

    Returns:
        Namespace: An object containing parsed arguments.
    '''

    parser = ArgumentParser()
    parser.add_argument('-i', '--input', type=Path, required=True, nargs='+',
                        help='Input PDB files to align.')
    parser.add_argument('-g', '--gui', action='store_true',
                        help='Run PyMOL instances with the GUI.')
    parser.add_argument('-r', '--ref', type=Path,
                        help='Reference PDB file to align to.')
    return parser.parse_args()

def align_structure():
    '''
    Align multiple PDB structures using PyMOL.

    Parses command-line arguments to get input files and alignment options,
    loads PDB files into PyMOL, aligns them to a reference structure (if provided),
    and saves the PyMOL session as a .pse file.

    Raises:
        FileNotFoundError: If any of the input PDB files or reference file cannot be found.
    '''

    # Handle script arguments
    user_args = script_args()

    # List of paths to pdb files
    pdb_files = user_args.input

    # Run the PyMOL instance with or without the GUI
    pymol_session = pymol.pymolPy3() if user_args.gui else pymol.pymolPy3(0)

    # Define the reference structure for alignment.
    # Defaults to the first pdb file from the inputs list
    ref_pdb_path = user_args.ref if user_args.ref else pdb_files[0]
    ref_pdb = ref_pdb_path.stem.replace(' ', '_')

    # Add reference structure to list of pdb files if it's not already
    if user_args.ref and user_args.ref not in pdb_files:
        pdb_files.append(user_args.ref)

    # Create a pymol command log file
    pymol_session(f'log_open {ref_pdb}_aligned.log')

    # Load each pdb file
    pymol_session(f'load {ref_pdb_path}')

    for pdb in pdb_files:
        pymol_session(f'load {pdb}')

    # Show all structures using cartoon representation
    pymol_session('hide all')
    pymol_session(f'set cartoon_transparency, 0.90,all')
    pymol_session('show cartoon,all')
    pymol_session(f'set cartoon_transparency, 0,{ref_pdb}')

    # Centre view on the reference structure
    pymol_session(f'center {ref_pdb}')

    # Align each structure to the reference structure
    for pdb in pdb_files:
        current_pdb = pdb.stem.replace(' ', '_')
        if current_pdb != ref_pdb:
            pymol_session(f'set cartoon_transparency, 0, {current_pdb}')
            # align is better for sequentially similar proteins
            # super or cealign is better for structurally similar proteins
            pymol_session(f'align {current_pdb}, {ref_pdb}')
            if user_args.gui:
                pymol_session('zoom all')
            pymol_session(f'set cartoon_transparency, 0.90, {current_pdb}')
        else:
            pass

    # Remove all water molecules
    pymol_session('remove resn hoh')

    # Center viewport view on the reference structure
    pymol_session(f'center {ref_pdb}')
    pymol_session(f'zoom {ref_pdb}')

    # Save current PyMOL session as a .pse file
    output_file = ref_pdb_path.parent / f'{ref_pdb}_align.pse'
    pymol_session(f'save {output_file}')
    print(f'Output file: {output_file}')

    for pdb in pdb_files:
        current_pdb = pdb.stem.replace(' ', '_')
        if current_pdb != ref_pdb:
            pymol_session(f'set cartoon_transparency, 0, {current_pdb}')
            pymol_session(f'cealign {ref_pdb}, {current_pdb}')
            if user_args.gui:
                pymol_session('zoom all')
            pymol_session(f'set cartoon_transparency, 0.90, {current_pdb}')
        else:
            pass

    pymol_session(f'center {ref_pdb}')
    pymol_session(f'zoom {ref_pdb}')

    output_file = ref_pdb_path.parent / f'{ref_pdb}_cealign.pse'
    pymol_session(f'save {output_file}')
    print(f'Output file: {output_file}')

    for pdb in pdb_files:
        current_pdb = pdb.stem.replace(' ', '_')
        if current_pdb != ref_pdb:
            pymol_session(f'set cartoon_transparency, 0, {current_pdb}')
            pymol_session(f'super {current_pdb}, {ref_pdb}')
            if user_args.gui:
                pymol_session('zoom all')
            pymol_session(f'set cartoon_transparency, 0.90, {current_pdb}')
        else:
            pass

    pymol_session(f'center {ref_pdb}')
    pymol_session(f'zoom {ref_pdb}')

    output_file = ref_pdb_path.parent / f'{ref_pdb}_super.pse'
    pymol_session(f'save {output_file}')
    print(f'Output file: {output_file}')


    pymol_session('log_close')

if __name__ == '__main__':
    align_structure()
