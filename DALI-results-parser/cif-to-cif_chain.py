from pathlib import Path
import warnings

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio import BiopythonWarning

import wget
from tqdm import tqdm

aligned_pdb_dir = Path('/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/Documents/GitHub/PhD-scripts/DALI-results-parser/results/1i8uA/pdb_dali_aligned')
pdb_dir = Path('/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/Documents/GitHub/PhD-scripts/DALI-results-parser/results/1i8uA/pdb_unaligned')

unique_pdbs = {}

for file_path in aligned_pdb_dir.iterdir():
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
    cif_file_path = pdb_dir / f'{cif}.cif'
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
            io = MMCIFIO()
            io.set_structure(chain_structure)
            single_chain_cif_path = pdb_dir / f'{cif}-{chain}.cif'
            io.save(str(single_chain_cif_path))
        except KeyError:
            # change to logging
            # print(f'Chain {chain} not found in {cif}.cif. Skipping...')
            continue
    