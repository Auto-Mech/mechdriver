""" drivers for coordinate scans
"""

def hrpotentials(tors_pots):
    """ Check hr pot to see if a new mimnimum is needed
    """

    print('\nHR potentials...')
    for name in tors_pots:

        print('- Rotor {}'.format(name))
        pot_str = ''
        for pot in tors_pots[name].values():
            pot_str += ' {0:.2f}'.format(pot)

        print('- Pot:{}'.format(pot_str))
