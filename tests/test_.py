""" Runs automech
"""

import os
import subprocess


PATH = os.getcwd()
INP_PATH = os.path.join(PATH, 'inp')
RUN_DAT_PATH = os.path.join(INP_PATH, 'run.dat')
CMD_LINE = (
    'python -u automech.py {0} >& {0}/run.log & disown %1'.format(PATH)
)


def test__automech():
    """ test automech.py
    """

    # Format the run.dat with the run-save dir paths
    with open(RUN_DAT_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(PATH))

    subprocess.call(CMD_LINE.split())


if __name__ == '__main__':
    test__automech()
