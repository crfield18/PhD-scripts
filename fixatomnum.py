'''
Fix the atom numbering for PDB files that have had new atoms added to them using
pdb4amber (AmberTools) https://ambermd.org/tutorials/basic/tutorial9/index.php, or
MolProbity http://molprobity.biochem.duke.edu.

Args:
    --input/-i: Path to one or more pdb files

Returns:
    Re-formatted PDB file following the naming convention "inputfilename"_output.pdb,
    located in the same directory as the input file.
'''

from pathlib import Path
from argparse import ArgumentParser

def script_args():
    '''
    Parse command-line arguments and return a Namespace object containing
    pathlib.Path paths to each user-defined input file.

    Returns:
        argparse.Namespace: A Namespace object containing input file paths.
    '''

    parser = ArgumentParser(description='Fix atom numbering in PDB file.')
    parser.add_argument('-i', '--input', type=Path, required=True, nargs='+',
                        help='Path to one or more pdb files')
    return parser.parse_args()

def line_fix_atom_num(line:str, atom_num_count:int):
    '''
    Format ATOM and TER lines from PDB files to fix the atom numbering in columns 7-11.

    Args:
        line:   String containing 1 line of an open PDB file.
        atom_num_count: Counter tracking the current atom number

    Returns:
        Tuple containing the updated line (str) and updated atom number count (int)
        (updated_line:str, updated_atom_num_count:int)
    '''

    if line.startswith('TER'):
        return 'TER', atom_num_count
    if line.startswith('ATOM'):
        return f'{line[0:6]}{str(atom_num_count).rjust(5)}{line[11:-1]}', atom_num_count + 1
    return line.strip('\n'), atom_num_count

def pdb_fix_atom_num(pdb_file):
    '''
    Formats each user-defined PDB file to fix the atom numbering in columns 7-11.

    Args:
        pdb_file:   pathlib.Path path to one of the input PDB files.

    Returns:
        None
    '''

    output_file_path =  pdb_file.parent / f'{pdb_file.stem}_output.pdb'

    with open(pdb_file, 'r', encoding='utf-8') as input_file, \
         open(output_file_path, 'w', encoding='utf-8') as output_file:
        atom_num_count = 1
        for line in input_file:
            newline, atom_num_count = line_fix_atom_num(line, atom_num_count)
            output_file.writelines(f'{newline}\n')

def main():
    '''
    Primary execution point of the script.

    Usage:
        To run this script, execute it from the command line as follows:
        $ python3 fixatomnum.py -i [path/to/file]
    '''

    user_args = script_args()
    for file in user_args.input:
        pdb_fix_atom_num(file)

if __name__ == '__main__':
    main()
