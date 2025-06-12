""" Runner
"""

from autorun._run import from_input_string as from_istring


SCRIPT_NAME = 'run_polyrate.sh'
INPUT_NAME = 'input.dat'
OUTPUT_NAMES = ('poly.fu6',)


def direct(script_str, run_dir, input_str, pot_str,
           script_name=SCRIPT_NAME,
           input_name=INPUT_NAME,
           output_names=OUTPUT_NAMES):
    """ Lazy runner polyrate
    """

    output_str = from_istring(
        script_str, run_dir, input_str,
        aux_dct={'input.fu40': pot_str},
        script_name=script_name,
        input_name=input_name,
        output_names=output_names)

    return input_str, output_str
