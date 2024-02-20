from pathlib import Path
import re
from argparse import ArgumentParser
import warnings

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio import BiopythonWarning

import wget
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
        results_index = self.dir_dali / 's001A_index.html'
        if not results_index.exists():
            wget.download(url=self.results_url, out=str(results_index), bar=None)

        with open(results_index, 'r', encoding='UTF8') as file:
            soup = BeautifulSoup(file, 'html.parser')
            self.pages = [a.get('href') for a in soup.find_all('a')]

        for page in self.pages:
            try:
                page_url = f'{self.results_url}{page}'
                output_file = self.dir_dali / page
                if not output_file.exists():

                    wget.download(url=page_url, out=str(output_file))
                else:
                    print(f'\r{output_file.name} found.\tSkipping...', end='', flush=True)
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
        for code, download_url in tqdm(self.dali_pdbs.items(),
                                       desc='',
                                       unit=' PDBs'):
            output_pdb = self.dir_pdb_aligned / f'{code}.pdb'
            if not output_pdb.exists():
                print(f'\n\rDownloading:\t{code}.pdb', end='', flush=True)
                wget.download(url=download_url, out=str(output_pdb))
                print(f'\n\rCleaning:\t{code}.pdb', end='', flush=True)
                html_to_pdb(output_pdb)
            else:
                print(f'\r{code}.pdb found. Skipping...', end='', flush=True)

    def cif_download(self):
        unique_pdbs = {}

        cif_originals_path = self.dir_pdb_clean / 'originals'
        cif_chains_path = self.dir_pdb_clean / 'chains'

        cif_originals_path.mkdir(parents=True, exist_ok=True)
        cif_chains_path.mkdir(parents=True, exist_ok=True)

        for file_path in self.dir_pdb_aligned.iterdir():
            if file_path.is_file():
                pdb_code = str(file_path.stem[0:4])
                chain_id = file_path.stem[4]

                if pdb_code not in unique_pdbs:
                    unique_pdbs[pdb_code] = list(chain_id)
                else:
                    unique_pdbs[pdb_code].append(chain_id)

        print(unique_pdbs)

        PDB_DL_URL_BASE = 'https://files.rcsb.org/download'

        warnings.filterwarnings('ignore', category=BiopythonWarning) # suppress PDBConstructionWarning
        parser = MMCIFParser()

        for cif, chain_list in tqdm(unique_pdbs.items(),
                                    desc='Downloading and Processing',
                                    unit=' mmCIF files'):
            cif_file_path = self.dir_pdb_clean / 'originals' / f'{cif}.cif'
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
                
                cif_chains_file_path = self.dir_pdb_clean / 'chains' / f'{cif}-chain.cif'
                if not cif_chains_file_path.exists():
                    try:
                        chain_structure = structure[0][chain]

                        io = MMCIFIO()
                        io.set_structure(chain_structure)
                        single_chain_cif_path = cif_chains_path / f'{cif}-{chain}.cif'
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
    dali_job.aligned_pdb_download()

    dali_job.cif_download()

if __name__ == '__main__':
    main()

