""" Runners for PAC99 program
"""

import os
import automol
import pac99_io.reader
from autorun._run import from_input_string


SCRIPT_NAME = 'run_pac99.sh'
INPUT_NAME = '{}.i97'
OUTPUT_NAMES = ('{}.o97', '{}.c97')

NEW_GROUPS_NAME = 'new.groups'


# Specialized runner
def nasa_polynomial(script_str, run_dir, input_str, name, formula,
                    convert=False):
    """ Generates NASA polynomial from run

        :param convert: convert the polynomial to more standard CHEMKIN
        :type convert
    """

    # Run PAC99 to get the output file
    formula_str = automol.form.string(formula)
    output_strs = direct(script_str, run_dir, input_str, formula_str)

    # Obtain the NASA polynomial, convert if necessary
    if output_strs is not None:
        c97_str = output_strs[1]
        poly_str = pac99_io.reader.nasa_polynomial(c97_str)
        if convert:
            poly_str = pac99_io.pac2ckin_poly(name, formula, poly_str)
    else:
        poly_str = None

    return poly_str


# Generalized runners
def direct(script_str, run_dir, input_str, formula_str):
    """ Generates an input file for a ThermP job runs it directly.

        Need formula input to run the script
        :param input_str: string of input file with .i97 suffix
    """

    aux_dct = {NEW_GROUPS_NAME: _new_groups_str()}

    input_name = INPUT_NAME.format(formula_str)
    output_names = tuple(name.format(formula_str) for name in OUTPUT_NAMES)
    output_strs = from_input_string(
        script_str, run_dir, input_str,
        aux_dct=aux_dct,
        script_name=SCRIPT_NAME,
        input_name=input_name,
        output_names=output_names)

    if not _check(output_strs):
        output_strs = None

    return output_strs


# Read the new groups file stored with src
def _new_groups_str():
    """ Read the new groups string
    """
    src_path = os.path.dirname(os.path.realpath(__file__))
    new_groups_path = os.path.join(src_path, 'aux', NEW_GROUPS_NAME)
    with open(new_groups_path, mode='r', encoding='utf-8') as fobj:
        new_groups_str = fobj.read()
    return new_groups_str


# Checks the output to see if there are errors
def _check(output_strs):
    """ assess the output (.o97, .c97 fileS)
    """

    o97_output_str, c97_output_str = output_strs

    success = True
    if 'INSUFFICIENT DATA' in o97_output_str:
        print('*ERROR: PAC99 fit failed, maybe increase temperature ranges?')
        success = False

    if not c97_output_str:
        print('No polynomial produced from PAC99 fits, check for errors')
        success = False

    return success
