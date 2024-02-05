# Create ray-traced images from Pymol Session (.pse) files

from pathlib import Path
from argparse import ArgumentParser
import pymolPy3 as pymol

def script_args():
    parser = ArgumentParser()
    parser.add_argument('-i', '--input',
                        type=Path, required=True,
                        nargs='+', help='Pymol session file(s) (.pse extension)')
    return parser.parse_args()

def main():
    pymol_session = pymol.pymolPy3(0)
    image_files = []
    input_pse = script_args().input
    for pse in input_pse:
        pymol_session('reinitialize')
        pymol_session(f'load {pse}')
        raw_figure = pse.with_suffix(".png")
        pymol_session(f'png {raw_figure}, width=3000, height=2000, dpi=300, ray=1')
        image_files.append(raw_figure)

if __name__ == '__main__':
    main()
