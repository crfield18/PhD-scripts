"""
Generates a template config file for use with AutoDock Vina/Vina-Carb
using a grid dimensions file generated with AutoDockTools.

Args:
    --input/-i: Path to one or more grid dimensions files

Returns:
    File named "config.txt" in the same directory as the input file
"""

from argparse import ArgumentParser
import pathlib

def script_args():
    """
    Parse command-line arguments and return a Namespace object containing
    pathlib.Path paths to each user-defined input file.

    Returns:
        argparse.Namespace: A Namespace object containing input file paths.
    """
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", type=pathlib.Path, required=True, nargs="+",
                        help="AutoDockTools gridbox output files")
    return parser.parse_args()

def parse_gridbox(gridbox_file):
    """
    Converts each user-defined input file into a template AutoDock Vina config file.

    Args:
        gridbox_file:   pathlib.Path path to one of the input files.

    Returns:
        None
    """

    gridbox_contents = []
    index_npts = None
    index_center = None

    with open(gridbox_file, "r", encoding="utf-8") as input_file:
        gridbox_contents = [line.split() for line in input_file if line.strip()]

    for index, item in enumerate(gridbox_contents):
        if item[0] == "npts":
            index_npts = index
        elif item[0] == "center":
            index_center = index

    variable_mappings = [
        ("center_x", 1),
        ("center_y", 2),
        ("center_z", 3),
        ("size_x", 1),
        ("size_y", 2),
        ("size_z", 3),
    ]

    with open(gridbox_file.with_name("config.txt"), "w", encoding="utf-8") as file:
        file.writelines("receptor = \nligand = \n\n")

        for variable, index in variable_mappings:
            if variable.startswith("center"):
                data_source = gridbox_contents[index_center]
            elif variable.startswith("size"):
                data_source = gridbox_contents[index_npts]

            file.writelines(f"{variable} = {data_source[index]}\n")

        file.writelines("\nenergy_range = 4\nexhaustiveness = 32\n")
        file.writelines("\nlog = \nout = \n")

def main():
    """
    Primary execution point of the script.

    Usage:
        To run this script, execute it from the command line as follows:
        $ python3 adv-config.py -i [path/to/file]
    """

    user_args = script_args()
    for file in user_args.input:
        parse_gridbox(file)

if __name__ == "__main__":
    main()
