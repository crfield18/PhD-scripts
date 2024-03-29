from pathlib import Path
import csv
from pymol import cmd
import psico.fullinit

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
                if method_function == cmd.cealign:
                    method_function(reference_object, obj)
                elif method_function in (cmd.align, cmd.super, cmd.tmalign):
                    method_function(obj, reference_object)
                else:
                    pass

                # Calculate and store RMSD and TM scores for the alignment
                rmsd = cmd.rms_cur(obj, reference_object)
                # tm = cmd.do(f'tmscore {obj}, {reference_object}')
                # tm = cmd.tmscore(obj, reference_object)
                # print(f'\r{obj}:\tRMSD:\t{rmsd:.5f} Å\tTM:{tm}    ', end='', flush=True)
                # print(f'\r{obj}:\tRMSD:\t{rmsd:.5f} Å    ', end='', flush=True)
                if obj not in results_dict:
                    results_dict[obj] = {}
                results_dict[obj]['RMSD'] = rmsd
                # results_dict[obj]['TM'] = tm
    # print(f'tm results type: {type(tm)}')
    print(results_dict)

    # # Write the alignment results of a csv file in the current working directory
    # with open(cwd / f'{reference_object}-{method}.csv', 'w', newline='') as csvfile:
    #     writer = csv.writer(csvfile)
    #     writer.writerow(['target', f'{method}_rmsd'])
    #     for key, value in results_dict.items():
    #         writer.writerow([key, value])

    # print(f'\nResults written to:\t{reference_object}-{method}.csv\n')

# Create callable function inside pymol called method_all for each alignment algorithm
def define_alignment_method(method_str):
    def align_all(reference_object):
        align_to_ref(reference_object, method=method_str)

    # function accepts 1 object to be used as the reference
    cmd.extend(f'{method_str}_all', align_all)
    cmd.auto_arg[0][f'{method_str}_all'] = [cmd.object_sc, 'object', '']

for algorithm in ['tmalign', 'align', 'cealign', 'super']:
    define_alignment_method(algorithm)


# Align using tmalign method
def tmalign_all(reference_object:str):
    align_to_ref(reference_object, method='tmalign')

# Align using align method
def align_all(reference_object:str):
    align_to_ref(reference_object, method='align')

# Align using cealign method
def cealign_all(reference_object:str):
    align_to_ref(reference_object, method='cealign')

# Align using super method
def super_all(reference_object:str):
    align_to_ref(reference_object, method='super')