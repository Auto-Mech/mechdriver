"""
 Generates NASA Polynomial from MESS+THERMP+PAC99 outputs
"""

import automol
import autorun
import ioformat
from mechlib.amech_io import writer
from mechlib.amech_io import printer as ioprinter


# NEW AUTORUN
def build_polynomial(spc_name, spc_dct, pf_path, nasa_path):
    """ Build a nasa polynomial
    """

    ioprinter.generating('NASA polynomials', nasa_path)

    # Generate forumula
    spc_dct_i = spc_dct[spc_name]
    formula_str = automol.inchi.formula_string(spc_dct_i['inchi'])
    formula_dct = automol.inchi.formula(spc_dct_i['inchi'])
    hform0 = spc_dct_i['Hfs'][0]

    thermp_script_str = autorun.SCRIPT_DCT['thermp']
    pac99_script_str = autorun.SCRIPT_DCT['pac99'].format(formula_str)

    # Copy MESSPF output file to THERMP run dir and rename to pf.dat
    pf_str = ioformat.pathtools.read_file(pf_path, 'pf.dat')

    hform298, poly_str = autorun.thermo.direct(
        thermp_script_str, pac99_script_str, nasa_path,
        pf_str, spc_name, formula_dct, hform0,
        enthalpyt=0.0, breakt=1000.0, convert=True)

    # Write the full CHEMKIN strings
    ckin_str = writer.ckin.nasa_polynomial(hform0, hform298, poly_str)
    full_ckin_str = '\n' + ckin_str

    return full_ckin_str
