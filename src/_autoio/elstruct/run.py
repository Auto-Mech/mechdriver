""" core run function
"""

from autorun import from_input_string


def direct(input_writer, script_str, run_dir, prog,
           geo, charge, mult, method, basis, **kwargs):
    """ Generates an input file for an electronic structure job and
        runs it directly.

        :param input_writer: elstruct writer module function for desired job
        :type input_writer: elstruct function
        :param script_str: string of bash script that contains
            execution instructions electronic structure job
        :type script_str: str
        :param run_dir: name of directory to run electronic structure job
        :type run_dir: str
        :param prog: electronic structure program to run
        :type prog: str
        :param geo: cartesian or z-matrix geometry
        :type geo: tuple
        :param charge: molecular charge
        :type charge: int
        :param mult: spin multiplicity
        :type mult: int
        :param method: electronic structure method
        :type method: str
        :returns: the input string, the output string, and the run directory
        :rtype: (str, str)
    """

    input_str = input_writer(
        prog=prog,
        geo=geo, charge=charge, mult=mult, method=method, basis=basis,
        **kwargs)

    output_strs = from_input_string(script_str, run_dir, input_str)
    output_str = output_strs[0]

    return input_str, output_str
