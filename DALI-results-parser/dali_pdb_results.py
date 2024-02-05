from pathlib import Path
import re
from argparse import ArgumentParser

import wget
import pymolPy3 as pymol
from bs4 import BeautifulSoup

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
    pymol_session = pymol.pymolPy3(0)

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
                    pdbid_with_chain = regex_finder(r'cd2=([a-zA-Z0-9]{5})', line)
                    pdb_url = regex_finder(r'http://[^"]+', line)

                    if pdbid_with_chain and pdb_url:
                        dali_pdbs[pdbid_with_chain.split('=')[-1]] = pdb_url

        self.dali_pdbs = dali_pdbs

        return self.dali_pdbs

    def pdb_download(self):
        for code, download_url in self.dali_pdbs.items():
            output_pdb = self.dir_pdb_aligned / f'{code}.pdb'
            if not output_pdb.exists():
                print(f'\rDownloading:\t{code}.pdb', end='', flush=True)
                wget.download(url=download_url, out=str(output_pdb))
                print(f'\rCleaning:\t{code}.pdb', end='', flush=True)
                html_to_pdb(output_pdb)
            else:
                print(f'{code}.pdb found. Skipping...')

def main():
    testrun = DaliResults()
    # Download all the DALI results pages
    testrun.download()
    # Extract all of the DALI-aligned PDB file download links
    testrun.pdb_links()
    # Download all the DALI-aligned PDB files
    # Maybe set an RMSD limit so you're not downloading ~5000 pdb files (e.g. 2.5 Å ?)
    # This should also limit memory usage form loading in pdb files to pymol
    testrun.pdb_download()

    # pymol_view(
    #     pdb_dir=testrun.dir_pdb_aligned,
    #     ref_pdb=Path('/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/Documents/GitHub/PhD-scripts/DALI/results/1i8uA/pdb_dali_aligned/1i8uA.pdb'),
    #     align_pdbs=False
    #     )

if __name__ == '__main__':
    main()
