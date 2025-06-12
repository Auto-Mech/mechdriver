""" tests varecof runners
"""

import os
import tempfile
# from ioformat import pathtools
import autorun

PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')
TMP_DIR = tempfile.mkdtemp()
print('dir', TMP_DIR)

# Structure/Coord information
TS_ZMA = ()
RCT_ZMAS = ()
BND_FRM_KEYS = ()
MEP_DISTANCES = ()

# Potential information
POTENTIALS = ()
NPOT = len(POTENTIALS)

# Other VRCTST info
VRC_DCT = {}

# Machine information
FORTRAN_COMPILER = 'gfortran'
MACHINE_DCT = {}

# Set the script strings
VARECOF_SCRIPT_STR = autorun.SCRIPT_DCT['varecof']
CONV_STRUCT_SCRIPT_STR = autorun.SCRIPT_DCT['varecof_conv_struct']
MCFLUX_SCRIPT_STR = autorun.SCRIPT_DCT['mcflux']

# Read the reference strings
# REF_FLUX_STR = ioformat.pathtools.read_file(DAT_PATH, 'flux.dat')


def test__():
    """ test autorun.varecof.write_varecof_input
    """

    # Generate the correction potentials
    autorun.varecof.compile_potentials(
        TMP_DIR, MEP_DISTANCES, POTENTIALS,
        BND_FRM_KEYS, FORTRAN_COMPILER,
        dist_restrict_idxs=(),
        pot_labels=(),
        pot_file_names=(),
        spc_name='mol')

    # Write the input strings
    inp_strs = autorun.varecof.write_input(
        TMP_DIR,
        TS_ZMA, RCT_ZMAS,
        NPOT, BND_FRM_KEYS,
        MACHINE_DCT, VRC_DCT)

    # Write the electronic structure input
    inp_strs += ()

    # Run VareCoF
    flux_str = autorun.varecof.flux_file(
        VARECOF_SCRIPT_STR, MCFLUX_SCRIPT_STR,
        TMP_DIR, inp_strs)
    print(flux_str)

    # Trim off the 0 energy (maybe do this in some function)

    # Check the fluxes to see if they are within some percentage threshold


if __name__ == '__main__':
    test__()
