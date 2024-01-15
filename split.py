'''
Determine the number of monovalent salt ions are required
for an explicit solvent MD simulation using the SPLIT method.
See: http://archive.ambermd.org/202002/0194.html

num_salt (No) = Expected number of salt ions
conc_salt (Co) = Salt concentration (M)
num_solv (Nw) = Number of water molecules
charge_solute (Q)  = Solute's charge

Returns:
    String in the form "addions2 pdb K+ {num_plus} Cl- {num_minus}"
    where num_plus = number of positively charged ions (K+) and
    num_minus = number of negatively charged ions (Cl-).

    For use with LEaP (AmberTools).
'''

from math import ceil

# User-editable values
input_co = 0.15
input_nw = 6817
input_q = 0

def is_even(num):
    '''
    Check whether the input number is even or odd.

    Args:
        num (int): Any integer.

    Returns:
        bool: True if the number is even, False if it's odd.
    '''
    if (num % 2) == 0:
        return True
    return False

def split_valid(num_salt, charge_solute):
    '''
    Check if the SPLIT method is valid for the provided data set.
    It is only valid when num_salt/charge_solute ≥ 1.

    Args:
        num_salt: Expected number of salt ions. (num_solv * (conc_salt / 56)
        charge_solute: Solute charge.

    Returns:
        bool: True if the system is valid, False otherwise.
    '''
    if charge_solute == 0 or num_salt/charge_solute >= 1:
        return True
    return False

def split(num_solv, conc_salt, charge_solute):
    '''
    Calculate the number of salt ions needed for a valid system.

    Args:
        num_solv (int): Number of solvent (water) molecules in the system.
        conc_salt (float): Salt concentration (M).
        charge_solute (float): Solute charge.

    Returns:
        tuple: A tuple containing the number of positively charged ions (K+) and
        the number of negatively charged ions (Cl-).
    '''
    num_salt = num_solv * (conc_salt / 56)
    if not split_valid(num_salt, charge_solute):
        return 'SPLIT not valid for this dataset (num_salt/charge_solute < 1)'

    num_plus = ceil(num_salt - (charge_solute / 2))
    num_minus = ceil(num_salt + (charge_solute / 2))

    return num_plus, num_minus

if __name__ == '__main__':
    num_k, num_cl = split(num_solv=input_nw, conc_salt=input_co, charge_solute=input_q)
    print(f'addions2 pdb K+ {num_k} Cl- {num_cl}')
