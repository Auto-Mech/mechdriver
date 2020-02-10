"""
Calculates Rate Constants from a Mechanism for given T,P
"""

import numpy as np
import chemkin_io


CKIN_FILES = [
#    './all_ls/l1.ckin',
#     './all_ls/l2.ckin',
     './all_ls/l3.ckin',
#     './all_ls/l3_cbs.ckin'
#    './abs_mig/l1_abs_mig.ckin',
#    './abs_mig/l2_abs_mig.ckin',
#    './abs_mig/l3_abs_mig.ckin',
#    './abs_mig/l3_cbs_abs_mig.ckin',
]

T_REF = 1.0
TEMPS = np.array([1000.0])
PRESSURES = [1.0, 10.0, 100.0]
PRINTPRES = 100.0


def calc_rates(mech1_str, t_ref, temps, pressures, typ=None):
    """ Read the reaction sections of a CHEMKIN files and
        calculate the rate constants
    """

    # Build dict for the rate constants
    mech1_reaction_block = chemkin_io.parser.util.clean_up_whitespace(
        chemkin_io.parser.mechanism.reaction_block(mech1_str))
    mech1_units = chemkin_io.parser.mechanism.reaction_units(mech1_str)
    mech1_ktp_dct = chemkin_io.calculator.rates.mechanism(
        mech1_reaction_block, mech1_units, t_ref, temps, pressures)

    # Get a list of the reactions aligining with the type
    reactions = _build_filtered_rxn_lst(mech1_ktp_dct, typ)
    reactions.sort(key=lambda x: x[0])

    # Print the reactions
    for reaction in reactions:
        if PRINTPRES in mech1_ktp_dct[reaction]:
            pressure = PRINTPRES
            print(_format_rxn_string(reaction), '      ', '{0:.2e}'.format(mech1_ktp_dct[reaction][pressure][0]))
        # else:
        #     pressure = 'high'
        # print(_format_rxn_string(reaction), '      ', '{0:.2e}'.format(mech1_ktp_dct[reaction][pressure][0]))
        # break


def total(mech1_str, t_ref, temps, pressures, typ=None):
    """ _ """
    # Build dict for the rate constants
    all_levels_ktp = []
    for ckin_file in CKIN_FILES:
        print(ckin_file)
        with open(ckin_file, 'r') as mech_file:
            mech_str = mech_file.read()
        mech_reaction_block = chemkin_io.parser.util.clean_up_whitespace(
            chemkin_io.parser.mechanism.reaction_block(mech_str))
        mech_units = chemkin_io.parser.mechanism.reaction_units(mech_str)
        mech_ktp_dct = chemkin_io.calculator.rates.mechanism(
            mech_reaction_block, mech_units, t_ref, temps, pressures)
        all_levels_ktp.append(mech_ktp_dct)

    # Get a list of the reactions aligining with the type
    all_levels_rxn = []
    for dct in all_levels_ktp:
        reactions = _build_filtered_rxn_lst(dct, typ)
        reactions.sort(key=lambda x: x[0])
        all_levels_rxn.append(reactions)

    # Build total dct of all levels
    total_dct = {}
    for i, rxn in enumerate(all_levels_rxn):
        for reaction in rxn:
            if PRINTPRES in reaction:
                pressure = PRINTPRES
            else:
                pressure = 'high'
            if reaction not in total_dct:
                total_dct[reaction] = _format_rxn_string(reaction)+'      '+'{0:.2e}'.format(all_levels_ktp[i][reaction][pressure][0])
            else:
                total_dct[reaction] += '      '+'{0:.2e}'.format(all_levels_ktp[i][reaction][pressure][0])

    for reaction, vals in total_dct.items():
        print(reaction)
        print(vals)



def calc_branching(ckin_files, t_ref, temps, pressures, typ=None):
    """ Read the reaction sections of CHEMKIN files and
        calculate the branching fractions
    """
    # Build dict for the rate constants
    mech_dct = {}
    for ckin_file in ckin_files:
        with open(ckin_file, 'r') as mech_file:
            mech_str = mech_file.read()
        reaction_block = chemkin_io.parser.util.clean_up_whitespace(
            chemkin_io.parser.mechanism.reaction_block(mech_str))
        units = chemkin_io.parser.mechanism.reaction_units(mech_str)
        kt_dct = chemkin_io.calculator.rates.mechanism(
            reaction_block, units, t_ref, temps, pressures)
        mech_dct.update(kt_dct)

    # Get a list of the reactions aligining with the type
    reactions = _build_filtered_rxn_lst(mech_dct, typ)

    # Get a dictionary of the total rate constants
    total_k_dct = {}
    for reaction in reactions:
        if PRINTPRES in mech_dct[reaction]:
            pressure = PRINTPRES
        else:
            pressure = 'high'
        if reaction[0] not in total_k_dct:
            total_k_dct[reaction[0]] = mech_dct[reaction][pressure][0]
        else:
            total_k_dct[reaction[0]] += mech_dct[reaction][pressure][0]

    # Loop over the mech dct to get the branching ratios
    print('mech dct')
    for reaction in reactions:
        if PRINTPRES in mech_dct[reaction]:
            pressure = PRINTPRES
        else:
            pressure = 'high'
        branch = (mech_dct[reaction][pressure][0] / total_k_dct[reaction[0]])
        if not np.isclose(branch, 1.0):
            print(reaction[0], reaction[1], '{0:.2f}'.format(float(branch)))
        # print(branch)


def _format_rxn_string(rxn):
    reacs, prods = rxn
    if len(prods) == 1:
        sep_sym = '<=>'
    else:
        sep_sym = '='
    reacs = ' + '.join(reacs)
    prods = ' + '.join(prods)
    rxn_str = '{0} {1} {2}'.format(reacs, sep_sym, prods)
    return rxn_str


def _build_filtered_rxn_lst(ktp_dct, typ):
    if typ == 'abstraction':
        nreacs, nprods = 2, 2
    elif typ == 'addition':
        nreacs, nprods = 2, 1
    elif typ == 'isomerization':
        nreacs, nprods = 1, 1

    reactions = []
    for rxn in ktp_dct:
        [reacs, prods] = rxn
        if len(reacs) == nreacs and len(prods) == nprods:
            reactions.append(rxn)

    return reactions


if __name__ == '__main__':
    print('rate constants')
    for ckin_file in CKIN_FILES:
        print(ckin_file)
        with open(ckin_file, 'r') as mech_file:
            mech_str = mech_file.read()
        calc_rates(mech_str, T_REF, TEMPS, PRESSURES, typ='abstraction')
    # print('\nbranching fraction')
    # calc_branching(CKIN_FILES, T_REF, TEMPS, PRESSURES, typ='abstraction')
    # total(CKIN_FILES, T_REF, TEMPS, PRESSURES, typ='abstraction')
