"""
compare thermo
"""

from chemkin_io.calculator import combine
from chemkin_io import plotter


def analyze_thermo(mech1_str, mech2_str, temps,
                   mech1_csv_str=None, mech2_csv_str=None,
                   dct_idx='name'):
    """ Read the thermo sections of two CHEMKIN files and
        produce plots of the thermochemical parameters together
    """

    # Build dictionaries containing the thermo data strings for each species
    if dct_idx == 'name':
        mech1_thermo_dct, mech2_thermo_dct = combine.build_thermo_name_dcts(
            mech1_str, mech2_str, temps)
    elif dct_idx == 'inchi':
        mech1_thermo_dct, mech2_thermo_dct = combine.build_thermo_inchi_dcts(
            mech1_str, mech2_str, mech1_csv_str, mech2_csv_str, temps)

    # Combine the two thermo dictionaries together under a common index
    thermo_vals_dct = combine.mechanism_thermo(
        mech1_thermo_dct, mech2_thermo_dct)

    # Build list of CHEMKIN mechanism names for plotting
    if dct_idx == 'name':
        names = [key for key in thermo_vals_dct]
    elif dct_idx == 'inchi':
        names = [combine.spc_name_from_inchi(mech1_csv_str, mech2_csv_str, key)
                 for key in thermo_vals_dct]

    # Plot the data from both mechanisms for each species
    plotter.thermo.build(thermo_vals_dct, temps, names)
