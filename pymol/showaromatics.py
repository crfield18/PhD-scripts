from pymol import cmd

def showaromatics(selection='all'):
    cmd.hide('everything', selection)
    cmd.show('cartoon', selection)
    cmd.color('gray', selection)
    cmd.select('aromatics', f'(resn tyr+trp+his) and {selection}')
    cmd.show('sticks', '(aromatics and !(name c+n+o))')
    cmd.color('magenta', 'aromatics')
    cmd.disable('aromatics')
    cmd.set('cartoon_smooth_loops', 0)

cmd.extend('showaromatics', showaromatics)
cmd.auto_arg[0]['showaromatics'] = [cmd.selection_sc, 'selection', '']
