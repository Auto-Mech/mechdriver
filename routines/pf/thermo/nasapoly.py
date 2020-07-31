"""
 Generates NASA Polynomial from MESS+THERMP+PAC99 outputs
"""

import os
import automol
import thermp_io
import pac99_io
from routines.pf import runner as pfrunner
from lib.amech_io import writer
from lib import pathtools


def build_polynomial(spc_name, spc_dct, temps,
                     pf_path, nasa_path, starting_path):
    """ Build a nasa polynomial
    """

    print('Generating NASA polynomials at path: {}'.format(nasa_path))

    # Generate forumula
    spc_dct_i = spc_dct[spc_name]
    formula = automol.inchi.formula_string(spc_dct_i['inchi'])
    formula_dct = automol.inchi.formula(spc_dct_i['inchi'])
    hform0 = spc_dct_i['Hfs'][0]

    # Go to NASA path
    if not os.path.exists(nasa_path):
        os.makedirs(nasa_path)
    pathtools.go_to(nasa_path)

    # Write and run ThermP to get the Hf298K and coefficients
    write_thermp_inp(formula, hform0, temps)
    pfrunner.run_thermp(pf_path, nasa_path)
    thermp_out_str = pathtools.read_file(nasa_path, 'thermp.out')
    hform298 = thermp_io.reader.hf298k(thermp_out_str)

    # Run PAC99 to get a NASA polynomial string in its format
    pfrunner.run_pac(formula, nasa_path)
    c97_file = pathtools.prepare_path(nasa_path, formula + '.c97')
    with open(c97_file, 'r') as file_obj:
        pac99_out_str = file_obj.read()
    # pac99_out_str = pathtools.read_file(o97_file)
    pac99_poly_str = pac99_io.reader.nasa_polynomial(pac99_out_str)

    # Obtain CHEMKIN string using PAC99 polynomial
    ckin_poly_str = pac99_io.pac2ckin_poly(
        spc_name, formula_dct, pac99_poly_str)

    # Write the full CHEMKIN strings
    header_str = '\n'
    nasa_str = writer.ckin.nasa_polynomial(hform0, hform298, ckin_poly_str)
    full_ckin_str = header_str + nasa_str

    print('\nCHEMKIN Polynomial:')
    print(full_ckin_str)

    # Go back to starting path
    pathtools.go_to(starting_path)

    return full_ckin_str


def write_thermp_inp(formula, hform0, temps,
                     enthalpyt=0.0, breakt=1000.0,
                     thermp_file_name='thermp.dat'):
    """ write the thermp input file
    """

    # Write thermp input file
    thermp_str = thermp_io.writer.input_file(
        ntemps=len(temps),
        formula=formula,
        delta_h=hform0,
        enthalpy_temp=enthalpyt,
        break_temp=breakt)

    # Write the file
    with open(thermp_file_name, 'w') as thermp_file:
        thermp_file.write(thermp_str)


def print_nasa_temps(temps):
    """ Print the polynomial fit temps
    """
    print('Attempting to fit NASA polynomials from',
          '200-1000 and 1000-3000 K ranges using\n',
          'temps from MESSPF file = {}.'.format(
              ' '.join(('{:.2f}'.format(x) for x in temps))))
