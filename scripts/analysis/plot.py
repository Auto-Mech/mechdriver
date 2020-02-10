"""
Create plots which compare thermochemical and kinetic values for two mechanisms
derived from the NASA polynomials and rate constants, respectively. Currently
assymed that the thermochemical and kinetic parameters will be read from
CHEMKIN-formatted mechanism files.
"""

import numpy as np
from chemkin_io.calculator import combine
from chemkin_io import plotter


CKIN1 = './natgas2/calcmech.txt'
CKIN2 = './natgas2/mechanism.txt'
CKIN1_SPC = './natgas2/calcspecies.csv'
CKIN2_SPC = './natgas2/species.csv'

with open(CKIN1, 'r') as ckin1:
    MECH1_STR = ckin1.read()
with open(CKIN2, 'r') as ckin2:
    MECH2_STR = ckin2.read()
with open(CKIN1_SPC, 'r') as ckin1:
    CSV1_STR = ckin1.read()
with open(CKIN2_SPC, 'r') as ckin2:
    CSV2_STR = ckin2.read()
T_REF = 1.0
TEMPS = np.arange(500.0, 2000.0, 25.00)
PRESSURES = [1.0, 10.0, 98.0]
RATEDIR = 'rates_natgas2'
THERMODIR = 'therm_plots'


def rates(mech1_str, mech2_str, t_ref, temps, pressures,
          mech1_csv_str=None, mech2_csv_str=None,
          dct_idx='name', plot_dir='rate_plots'):
    """ Read the reaction sections of two CHEMKIN files and
        produce plots of the rate constants together
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
        mech1_ktp_dct, mech2_ktp_dct, temps,
        mech2_thermo_dct=mech2_thermo_dct)

    # Build list of CHEMKIN mechanism names for plotting
    if dct_idx != 'name':
        ktp_dct = combine.conv_ich_to_name_ktp_dct(ktp_dct, mech1_csv_str)
    names = [key for key in ktp_dct]

    # Plot the data from both mechanisms for each species
    plotter.rates.build(
        ktp_dct, temps, plot_dir=plot_dir, names=names)


def thermo(mech1_str, mech2_str, temps,
           mech1_csv_str=None, mech2_csv_str=None,
           dct_idx='name', plot_dir='therm_plots'):
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
    plotter.thermo.build(
        thermo_vals_dct, temps, plot_dir=plot_dir, names=names)


if __name__ == '__main__':
    rates(MECH1_STR, MECH2_STR, T_REF, TEMPS, PRESSURES,
          mech1_csv_str=CSV1_STR, mech2_csv_str=CSV2_STR,
          dct_idx='inchi', plot_dir=RATEDIR)
