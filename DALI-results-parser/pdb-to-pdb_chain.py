from pathlib import Path
import wget

pdb_dir = Path('/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/Documents/GitHub/PhD-scripts/DALI/results/1i8uA/pdb_dali_aligned')

unique_pdbs = {}

for file_path in pdb_dir.iterdir():
    if file_path.is_file():
        pdb_code = str(file_path.stem[0:4])
        chain_id = file_path.stem[4]

        print(pdb_code, chain_id)
        if pdb_code not in unique_pdbs:
            unique_pdbs[pdb_code] = list(chain_id)
        else:
            unique_pdbs[pdb_code].append(chain_id) 

pdb_dl_url_base = 'https://files.rcsb.org/download'

for pdb, chain_list in unique_pdbs:
    wget.download(url=f'{pdb_dl_url_base}/{pdb}.pdb', out=str(pdb_dir / f'{pdb}.pdb'))
    for chain in chain_list:
        with open(pdb_dir / f'{pdb}.pdb', 'r') as full_pdb, open(pdb_dir / f'{pdb}-{chain}.pdb') as single_chain_pdb:
            for line in full_pdb:
                pass



class PDBInterpret:
    def __init__(self) -> None:
        self.line = ''
        self.record = self.line[0:6]
        self.res_name = self.line[17:20]
        self.chain_id = self.line[21]
    
    def line_filter(self, chain_to_keep:str):
        if self.record in ('ATOM  ', 'HETATM') and self.chain_id != chain_to_keep:
            return None
        if self.res_name == 'HOH':
            return None
        return self.line
