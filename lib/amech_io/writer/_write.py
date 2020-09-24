""" Write a series of input for some program
"""

import os
import stat


def write_files(inp_files, run_path, exe_names=()):
    """ Write all of the VaReCoF inut files and run the code
    """

    # Write all of the VaReCoF input files
    for (inp_str, inp_name) in inp_files:
        file_name = os.path.join(run_path, inp_name)
        print(file_name)
        # Write the file
        with open(file_name, 'w') as inp_file:
            inp_file.write(inp_str)
        # Make file an executable
        if inp_name in exe_names:
            os.chmod(file_name, mode=os.stat(file_name).st_mode | stat.S_IEXEC)
