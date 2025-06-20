""" Tests the writing of the energy transfer section
"""

import os
from ioformat import pathtools
import onedmin_io.writer


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')


# Input parameters
RANSEED = 28384920
NSAMP = 10
TARGET_XYZ_NAME = 'target.xyz'
BATH_XYZ_NAME = 'bath.xyz'
SMIN = 3.77945
SMAX = 9.44863

# Submission script parameters
NJOBS = 20
RUN_DIR = 'job'
EXE_PATH = 'onedmin.x'


def test__input_file():
    """ test onedmin_io.writer.input_file
    """

    onedmin_inp1_str = onedmin_io.writer.input_file(
        NSAMP, SMIN, SMAX, ranseed=RANSEED)

    onedmin_inp2_str = onedmin_io.writer.input_file(
        NSAMP, SMIN, SMAX, ranseed=RANSEED,
        target_xyz_name=TARGET_XYZ_NAME, bath_xyz_name=BATH_XYZ_NAME,
        spin_method=1)

    assert onedmin_inp1_str == pathtools.read_file(DAT_PATH, 'onedmin1.inp')
    assert onedmin_inp2_str == pathtools.read_file(DAT_PATH, 'onedmin2.inp')

    # Test for random seed
    onedmin_inp3_str = onedmin_io.writer.input_file(
        NSAMP, SMIN, SMAX, ranseed=None)
    lines = onedmin_inp3_str.strip().splitlines()

    # Convert first line to int (should work if ranseed is correct)
    # Check if it falls in allowed range set in function
    seed = int(lines[0])
    assert 1.0e7 < seed < 9.9e7

    # Check remaining lines
    assert lines[1] == '10'
    assert lines[2] == 'target.xyz'
    assert lines[3] == 'bath.xyz'
    assert lines[4] == '2.000 5.000'
    assert lines[5] == '2'


def test__submission_script():
    """ test onedmin_io.writer.submission_script
    """

    subscript_str = onedmin_io.writer.submission_script(
        NJOBS, RUN_DIR, EXE_PATH)

    assert subscript_str == pathtools.read_file(DAT_PATH, 'subscript.inp')
