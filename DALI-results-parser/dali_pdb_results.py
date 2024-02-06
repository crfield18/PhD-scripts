from pathlib import Path
import re
from argparse import ArgumentParser
import warnings

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio import BiopythonWarning

import wget
import pymolPy3 as pymol
from bs4 import BeautifulSoup
from tqdm import tqdm

def script_args():
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, help='DALI results URL')
    return parser.parse_args()

def regex_finder(regex_pattern, url: str):
    match = re.search(regex_pattern, url)
    if match:
        return match.group()
    return None

def html_pdb_line_filter(line: str):
    if not line.strip() or '<' in line:
        return None
    return line

def html_to_pdb(pdb_file):
    with open(pdb_file, 'r+', encoding='UTF8') as file:
        lines = file.readlines()
        # Remove any non-PDB format lines from the file
        filtered_lines = filter(lambda line: html_pdb_line_filter(line), lines)
        # Move the cursor to the beginning of the file
        file.seek(0)
        file.writelines(filtered_lines)
        file.truncate()

def pymol_view(pdb_dir:Path, ref_pdb:Path, align_pdbs:bool):
    pymol_session = pymol.pymolPy3()

    loaded_pdbs = []

    for file_path in pdb_dir.iterdir():
        if file_path.is_file():
            # Perform operations with the file here
            pymol_session(f'load {file_path}')
            loaded_pdbs.append(file_path.stem)

    # Show all structures using cartoon representation
    pymol_session('hide all')
    pymol_session('show cartoon,all')

    # Centre view on the reference structure
    pymol_session(f'center {ref_pdb}')
    pymol_session('zoom all')

    # Align PDBs to reference
    if align_pdbs:
        output_file = f'{ref_pdb}_aligned-pymol.pse'
        for pdb in loaded_pdbs:
            pymol_session(f'align {pdb}, {ref_pdb}')
    else:
        output_file = f'{ref_pdb}_aligned-dali.pse'

    # Save Pymol session
    pymol_session(f'save {output_file}')
    print(f'Output file: {output_file}')

class DaliResults:
    def __init__(self) -> None:
        self.results_url = script_args().input
        self.job_id = self.results_url.split('/')[-2]

        # Set up download directory structure
        self.cwd = Path.cwd()
        self.dir_downloads = self.cwd / 'results' / self.job_id

        directories = {
            'dir_dali': 'dali_results_pages',
            'dir_pdb_aligned': 'pdb_dali_aligned',
            'dir_pdb_clean': 'pdb_unaligned'
        }

        for name, dirname in directories.items():
            setattr(self, name, self.dir_downloads / dirname)
            getattr(self, name).mkdir(parents=True, exist_ok=True)

    def download(self):
        results_index = self.dir_dali / f'{self.job_id}_index.html'
        if not results_index.exists():
            wget.download(url=self.results_url, out=str(results_index))

        with open(results_index, 'r', encoding='UTF8') as file:
            soup = BeautifulSoup(file, 'html.parser')
            self.pages = [a.get('href') for a in soup.find_all('a')]

        for page in self.pages:
            try:
                page_url = f'{self.results_url}{page}'
                output_file = self.dir_dali / page
                if not output_file.exists():
                    print(f'Downloading {output_file}')
                    wget.download(url=page_url, out=str(output_file))
                else:
                    print(f'{output_file.name} found.\tSkipping...')
            except Exception as e:
                print(f'\nError downloading {output_file.name}: {e}.')

        return self.pages

    def pdb_links(self):
        dali_pdbs = {}
        with open(self.dir_dali / f'{self.job_id}.html', 'r', encoding='UTF8') as results_html:
            for line in results_html:
                if 'VALUE="cd2=' not in line:
                    pass
                else:
                    # Extract the PDB code + chain and the URL to the aligned pdb file from results html
                    pdbid_with_chain = regex_finder(r'cd2=([0-9][a-zA-Z0-9]{4})', line)
                    pdb_url = regex_finder(r'http://[^"]+', line)

                    if pdbid_with_chain and pdb_url:
                        dali_pdbs[pdbid_with_chain.split('=')[-1]] = pdb_url

        self.dali_pdbs = dali_pdbs

        return self.dali_pdbs

    def aligned_pdb_download(self):
        for code, download_url in tqdm(self.dali_pdbs.items(), desc='Downloading DALI aligned PDBs', unit=' PDBs'):
            output_pdb = self.dir_pdb_aligned / f'{code}.pdb'
            if not output_pdb.exists():
                print(f'\rDownloading:\t{code}.pdb', end='', flush=True)
                wget.download(url=download_url, out=str(output_pdb))
                print(f'\rCleaning:\t{code}.pdb', end='', flush=True)
                html_to_pdb(output_pdb)
            else:
                print(f'{code}.pdb found. Skipping...')

    def cif_download(self):
        unique_pdbs = {}

        for file_path in self.pdb_dali_aligned.iterdir():
            if file_path.is_file():
                pdb_code = str(file_path.stem[0:4])
                chain_id = file_path.stem[4]

                if pdb_code not in unique_pdbs:
                    unique_pdbs[pdb_code] = list(chain_id)
                else:
                    unique_pdbs[pdb_code].append(chain_id)

            # print(unique_pdbs)

            PDB_DL_URL_BASE = 'https://files.rcsb.org/download'

            warnings.filterwarnings('ignore', category=BiopythonWarning) # suppress PDBConstructionWarning
            parser = MMCIFParser()

            for cif, chain_list in tqdm(unique_pdbs.items(),
                                        desc='Downloading and Processing',
                                        unit=' mmCIF files'):
                cif_file_path = self.pdb_unaligned / f'{cif}.cif'
                if not cif_file_path.exists():
                    try:
                        # Download the PDB file
                        wget.download(url=f'{PDB_DL_URL_BASE}/{cif}.cif', out=cif_file_path, bar=None)
                    except Exception as e:
                        # print(f'\nError downloading {cif}: {e}.\n')
                        # change this print to some kind of logging thing??
                        continue

                structure = parser.get_structure(structure_id=cif, filename=cif_file_path)

                for chain in tqdm(chain_list,
                                desc=f'Processing Chains for {cif}',
                                unit=' chains',
                                leave=False):
                    try:
                        chain_structure = structure[0][chain]

                        # # Cut HETATM entries
                        # chain_structure = chain_structure.copy()
                        # heteroatoms = [residue for residue in chain_structure.get_residues() if residue.id[0] != " "]
                        # for heteroatom in heteroatoms:
                        #     chain_structure.detach_child(heteroatom.id)

                        # # Skip chain if all HETATM
                        # if not chain_structure.child_list:
                        #     continue

                        io = MMCIFIO()
                        io.set_structure(chain_structure)
                        single_chain_cif_path = self.pdb_unaligned / f'{cif}-{chain}.cif'
                        io.save(str(single_chain_cif_path))
                    except KeyError:
                        # change to logging
                        # print(f'Chain {chain} not found in {cif}.cif. Skipping...')
                        continue

def main():
    dali_job = DaliResults()
    # Download all the DALI results pages
    dali_job.download()
    # Extract all of the DALI-aligned PDB file download links
    dali_job.pdb_links()
    # Download all the DALI-aligned PDB files
    # Maybe set an RMSD limit so you're not downloading ~5000 pdb files (e.g. 2.5 Å ?)
    # This should also limit memory usage form loading in pdb files to pymol
    dali_job.aligned_pdb_download()

    # # Load all DALI-aligned PDB files into a Pymol session
    # pymol_view(
    #     pdb_dir=dali_job.dir_pdb_aligned,
    #     ref_pdb=Path('/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/Documents/GitHub/PhD-scripts/DALI/results/1i8uA/pdb_dali_aligned/1i8uA.pdb'),
    #     align_pdbs=False
    #     )

    dali_job.cif_download()

    # # Load all clean PDB files (.mmcif format to deal with large models) into a new Pymol session and align
    # # Need to add method for capturing RMSD values from pymol to construct a df with DALI and pymol RMSD values
    # pymol_view(
    # pdb_dir=Path('/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/Documents/GitHub/PhD-scripts/DALI-results-parser/results/1i8uA/pdb_unaligned/chains'),
    # ref_pdb=Path('/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/Documents/GitHub/PhD-scripts/DALI-results-parser/results/1i8uA/pdb_unaligned/chains/1i8u-A.cif'),
    # align_pdbs=True
    # )

if __name__ == '__main__':
    main()
