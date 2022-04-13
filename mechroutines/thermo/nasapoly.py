"""
 Generates NASA Polynomial from MESS+THERMP+PAC99 outputs
"""

import automol
import autorun
import ioformat
from mechlib.amech_io import writer
from mechlib.amech_io import printer as ioprinter


def build_polynomial(spc_name, spc_dct, pf_path, nasa_path,
                     spc_locs_idx=None, spc_mod=None):
    """ For a given species: obtain partition function data read from a
        MESSPF output file currently existing in the RUN filesystem as well as
        the previously computed 0 K heat-of-formation from the species
        dictionary. Then use this to run ThermP and PAC99 in the RUN filesystem
        to generate a ChemKin-formatted 7-coefficient NASA polynomial.

        :param spc_name: mechanism name of species to write MESSPF input for
        :type spc_name: str
        :param spc_dct:
        :type spc_dct:
        :param pf_path: path to existing MESSPF file in RUN filesystem
        :type pf_path: str
        :param nasa_path: path to run ThermP+PAC99 in RUN filesystem
        :type nasa_path: str
        :rtype: str
    """

    # Read the temperatures from the pf.dat file, check if viable
    ioprinter.nasa('fit', path=pf_path)
    ioprinter.generating('NASA polynomials', nasa_path)

    # Generate forumula
    spc_dct_i = spc_dct[spc_name]
    formula_str = automol.chi.formula_string(spc_dct_i['inchi'])
    formula_dct = automol.chi.formula(spc_dct_i['inchi'])

    if spc_locs_idx == 'final':
        hform0 = spc_dct_i['Hfs']['final'][0]
        spc_label = spc_name
    elif spc_locs_idx is not None:
        hform0 = spc_dct_i['Hfs'][spc_locs_idx][spc_mod][0]
        spc_label = spc_name  # + '_{:g}'.format(spc_locs_idx)
        # spc_label = spc_name + '_' + spc_locs[1][:5]
    else:
        hform0 = spc_dct_i['Hfs'][0]
        spc_label = spc_name
    thermp_script_str = autorun.SCRIPT_DCT['thermp']
    pac99_script_str = autorun.SCRIPT_DCT['pac99'].format(formula_str)

    # Copy MESSPF output file to THERMP run dir and rename to pf.dat
    pf_str = ioformat.pathtools.read_file(pf_path, 'pf.dat')
    hform298, poly_str = autorun.thermo(
        thermp_script_str, pac99_script_str, nasa_path,
        pf_str, spc_label, formula_dct, hform0,
        enthalpyt=0.0, breakt=1000.0, convert=True)

    # Write the full CHEMKIN strings
    ckin_str = '\n' + writer.ckin.nasa_polynomial(hform0, hform298, poly_str)

    return ckin_str
