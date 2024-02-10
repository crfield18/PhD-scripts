from pathlib import Path
import csv
from pymol import cmd

cwd = Path.cwd()

# Get a list of every object present in the current session (excluding 'sele')
def non_sele_objs():
    all_objects = cmd.get_object_list()
    non_sele_obj_list = [obj for obj in all_objects if obj != 'sele']
    return non_sele_obj_list

# Align each object (excluding 'sele') to reference_object using the specified method
def align_to_ref(reference_object:str, method:str):
    # Centre viewport on reference object
    cmd.center(reference_object)

    # Assign method strings to pymol commands
    alignment_functions = {
        'align': cmd.align,
        'cealign': cmd.cealign,
        'super': cmd.super,
        # Requires TMalign executable. See: https://pymolwiki.org/index.php/TMalign
        'tmalign': cmd.tmalign
    }

    results_dict = {}

    # Run command for each object (excluding 'sele' and the chosen reference object)
    method_function = alignment_functions.get(method)

    if method_function:
        print(f'\nReference object:\t{reference_object}\n')
        for obj in non_sele_objs():
            if obj != reference_object:
                # Extract the RMSD value from the different result types
                results = method_function(obj, reference_object)
                print(type(results))
                print(results)
                if method in ('super', 'align'):
                    rmsd = results[0]
                elif method == 'cealign':
                    rmsd = results["RMSD"]
                print(f'\r{obj}:\t{rmsd:.5f} Å    ', end='', flush=True)
                results_dict[obj] = rmsd

    # Write the alignment results of a csv file in the current working directory
    with open(cwd / f'{reference_object}-{method}.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['target', f'{method}_rmsd'])
        for key, value in results_dict.items():
            writer.writerow([key, value])

    print(f'\nResults written to:\t{reference_object}-{method}.csv\n')

# Create callable function inside pymol called method_all for each alignment algorithm
def define_alignment_method(method_str):
    def align_all(reference_object):
        align_to_ref(reference_object, method=method_str)

    # function accepts 1 object to be used as the reference
    cmd.extend(f'{method_str}_all', align_all)
    cmd.auto_arg[0][f'{method_str}_all'] = [cmd.object_sc, 'object', '']

for algorithm in ['tmalign', 'align', 'cealign', 'super']:
    define_alignment_method(algorithm)
