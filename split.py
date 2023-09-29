'''Determine the number of monovalent salt ions are required
for an explicit solvent MD simulation using the SPLIT method
(http://archive.ambermd.org/202002/0194.html)

No = Expected number of salt ions
Co = Salt concentration (M)
Nw = Number of water molecules
Q  = Solute's charge'''

from math import ceil

# User editable values
Input_Co = 0.15
Input_Nw = 6817
Input_Q = 0

def is_even(num:int):
    '''Checks whether the input num is even or odd.'''
    if (num % 2) == 0:
        return True
    return False

def split_valid(No:float, Q:float):
    '''Checks if the SPLIT method is valid for the provided data set.
    It is only valid when No/Q ≥ 1'''
    if Q == 0 or No/Q >= 1:
        return True
    return False

def split(Nw:float, Co:float, Q:float):
    '''Calculates the number of salt ions needed'''
    # if is_even(Q) == False:
    #     Q += 1
    No = Nw*(Co/56)
    if not split_valid(No, Q):
        return 'SPLIT not valid for this dataset (No/Q < 1)'

    nplus = ceil(No - (Q/2))
    nminus = ceil(No + (Q/2))

    return nplus, nminus

if __name__ == '__main__':
    nplus, nminus = split(Nw=Input_Nw, Co=Input_Co, Q=Input_Q)
    print(f'addions2 pdb K+ {nplus} Cl- {nminus}')
