""" Runs automech
"""

import os
import shutil
import tempfile
import subprocess


# Set path where test input file exist
PATH = os.path.dirname(os.path.realpath(__file__))
INP_DIR = os.path.join(PATH, 'inp')

# Set paths where tests will run
TMP_DIR = tempfile.mkdtemp()
RUN_DIR = os.path.join(TMP_DIR, 'run')
SAVE_DIR = os.path.join(TMP_DIR, 'save')
print('dir', TMP_DIR)

# Set names of input
INP_NAMES = ('run.template', 'mechanism.dat', 'species.csv',
             'models.dat', 'theory.dat')

# Set command line
EXE_PATH = os.path.join(PATH, '../bin/automech.py')
CMD_LINE = 'python -u {0} {1} & disown %1'.format(EXE_PATH, TMP_DIR)


def test__():
    """ test automech.py
    """
    _copy_input()
    _rewrite_run_dat()

    logfile = open('{0}/run.log'.format(TMP_DIR), 'w')
    subprocess.call(CMD_LINE.split(), stdout=logfile, stderr=logfile)


def _copy_input():
    """ Copy all of the input files
    """

    inp_paths = tuple(
        os.path.join(INP_DIR, name) for name in INP_NAMES)
    tmp_paths = tuple(
        os.path.join(TMP_DIR, 'inp', name) for name in INP_NAMES)
    os.makedirs(
        os.path.join(TMP_DIR, 'inp'))

    for inp_path, tmp_path in zip(inp_paths, tmp_paths):
        shutil.copy(inp_path, tmp_path)


def _rewrite_run_dat():
    """ Read the run.dat and replace run_prefix and save_prefix
    """
    # Read the input file into a string
    inp_file = os.path.join(TMP_DIR, 'inp', 'run.template')
    with open(inp_file, 'r') as fobj:
        inp_str = fobj.read()

    # Put the dynamically loaded tmp paths into the input
    new_inp_str = inp_str.format(RUN_DIR, SAVE_DIR)

    # Write the input
    new_inp_file = os.path.join(TMP_DIR, 'inp', 'run.dat')
    with open(new_inp_file, 'w') as fobj:
        inp_str = fobj.write(new_inp_str)


if __name__ == '__main__':
    test__()
