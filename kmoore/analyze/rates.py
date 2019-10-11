"""
compare thermo
"""

from chemkin_io.calculator import combine
from chemkin_io import plotter


def analyze_rates(mech1_str, mech2_str, t_ref, temps, pressures,
                  mech1_csv_str=None, mech2_csv_str=None,
                  dct_idx='name'):
    """ Read the reaction sections of two CHEMKIN files and
        produce plots of the thermochemical parameters together
    """

    # Build dictionaries containing the thermo data strings for each species
    if dct_idx == 'name':
        _, mech2_thermo_dct = combine.build_thermo_name_dcts(
            mech1_str, mech2_str, temps)
        mech1_ktp_dct, mech2_ktp_dct = combine.build_reaction_name_dcts(
            mech1_str, mech2_str,
            t_ref, temps, pressures)
    elif dct_idx == 'inchi':
        _, mech2_thermo_dct = combine.build_thermo_inchi_dcts(
            mech1_str, mech2_str, mech1_csv_str, mech2_csv_str, temps)
        mech1_ktp_dct, mech2_ktp_dct = combine.build_reaction_inchi_dcts(
            mech1_str, mech2_str, mech1_csv_str, mech2_csv_str,
            t_ref, temps, pressures)


    # Combine the two thermo dictionaries together under a common index
    ktp_dct = combine.mechanism_rates(
        mech1_ktp_dct, mech2_ktp_dct,
        mech2_thermo_dct,
        temps)

    # Build list of CHEMKIN mechanism names for plotting
    if dct_idx == 'name':
        names = [key for key in ktp_dct]
    # elif dct_idx == 'inchi':
    #     names = [combine.rxn_name_from_inchi(mech1_csv_str, mech2_csv_str, key)
    #              for key in ktp_dct]

    # Plot the data from both mechanisms for each species
    plotter.rates.build(ktp_dct, temps, names)
