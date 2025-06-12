""" Runner
"""

import thermp_io.writer
from autorun._run import from_input_string


SCRIPT_NAME = 'run_thermp.sh'
INPUT_NAME = 'thermp.dat'
OUTPUT_NAMES = ['thermp.out', '{}.i97']


# Generalized runner
def direct(script_str, run_dir,
           pf_str, formula_str, hform0, temps,
           enthalpyt=0.0, breakt=1000.0):
    """ Generates an input file for a ThermP job runs it directly.
    """

    input_str = thermp_io.writer.input_file(
        ntemps=len(temps),
        formula=formula_str,
        delta_h=hform0,
        enthalpy_temp=enthalpyt,
        break_temp=breakt)
    aux_dct = {'pf.dat': pf_str}

    output_names = OUTPUT_NAMES.copy()
    output_names[1] = output_names[1].format(formula_str)
    output_strs = from_input_string(
        script_str, run_dir, input_str,
        aux_dct=aux_dct,
        script_name=SCRIPT_NAME,
        input_name=INPUT_NAME,
        output_names=output_names)

    return input_str, output_strs
