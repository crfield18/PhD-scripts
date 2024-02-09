from pymol import cmd

def non_sele_objs():
    # Get a list of every object present in the current session (excluding 'sele')
    all_objects = cmd.get_object_list()
    non_sele_obj_list = [obj for obj in all_objects if obj != 'sele']
    return non_sele_obj_list

def align_to_ref(reference_object:str, method:str):
    # Centre viewport on reference object
    cmd.center(reference_object)

    # Assign method strings to pymol commands
    alignment_functions = {
        'align': cmd.align,
        'cealign': cmd.cealign,
        'super': cmd.super
        # # Requires TMalign executable. See: https://pymolwiki.org/index.php/TMalign
        # 'tmalign': cmd.tmalign
    }

    # Run command for each object (excluding 'sele' and the chosen reference object)
    method_function = alignment_functions.get(method)
    if method_function:
        print(f'\nReference object:\t{reference_object}\n')
        for obj in non_sele_objs():
            if obj != reference_object:
                results = method_function(obj, reference_object)
                # super and align
                if type(results) == tuple:
                    rmsd = results[0]
                # cealign
                elif type(results) == dict:
                    rmsd = results["RMSD"]
                print(f'{obj}:\t{rmsd:.5f} Å')
    else:
        pass

# Create callable function inside pymol called method_all for each alignment algorithm
def define_alignment_method(method_str):
    def align_all(reference_object):
        align_to_ref(reference_object, method=method_str)
        print(f'{method_str}_all complete!')
    # function accepts 1 object to be used as the reference
    cmd.extend(f'{method_str}_all', align_all)
    cmd.auto_arg[0][f'{method_str}_all'] = [cmd.object_sc, 'object', '']

alignment_methods = ['align', 'cealign', 'super']

for algorithm in alignment_methods:
    define_alignment_method(algorithm)
