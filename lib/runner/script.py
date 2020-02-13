""" Run bash scripts
"""

import os
import subprocess
import warnings
import stat


def run_script(script_str, run_dir):
    """ run a program from a script
    """

    script_name = 'build.sh'
    with _EnterDirectory(run_dir):
        # Write the submit script to the run directory
        try:
            os.remove('build.sh')
            with open(script_name, 'w') as script_obj:
                script_obj.write(script_str)
        except IOError:
            with open(script_name, 'w') as script_obj:
                script_obj.write(script_str)

        # Make the script executable
        os.chmod(script_name, mode=os.stat(script_name).st_mode | stat.S_IEXEC)

        # Call the program
        try:
            subprocess.check_call('./{:s}'.format(script_name))
        except subprocess.CalledProcessError:
            # If the program failed, continue with a warning
            warnings.warn("run failed in {}".format(run_dir))


class _EnterDirectory():

    def __init__(self, directory):
        assert os.path.isdir(directory)
        self.directory = directory
        self.working_directory = os.getcwd()

    def __enter__(self):
        os.chdir(self.directory)

    def __exit__(self, _exc_type, _exc_value, _traceback):
        os.chdir(self.working_directory)


# LIBRARY OF STANDARD INPUT STRINGS
def read_inp_scripts(script_dir):
    """ Read in the scripts the user submits as a strings
    """

    # Build a dct for the scripts string
    script_dct = {}
    script_names = [name for name in os.listdir(script_dir)]
    for name in script_names:
        with open(os.path.join(scripts_dir, name), 'r') as script_file:
            script_dct[name.replace('.sh', '')] = script_file.read()

    # Check the submission strings for viability
    # assert check_script_dct 
    return script_dct


# def std_script(prog):
#     """ Set the name of a standard script
#     """
#     if proj == 'projrot':
#         sub_str = PROJROT
#     return sub_str

# ProjRot
PROJROT = ("#!/usr/bin/env bash\n"
           "/lcrc/project/CMRP/pacc/ProjRot/build/RPHt.exe >& /dev/null")

# MESS
MESSPF = ("#!/usr/bin/env bash\n"
          "export OMP_NUM_THREADS=10\n"
          "messpf pf.inp pf.out >> stdout.log &> stderr.log")
MESSRATE = ("#!/usr/bin/env bash\n"
            "export OMP_NUM_THREADS=10\n"
            "mess mess.inp rate.out >> stdout.log &> stderr.log")

# VaReCoF
VARECOF = ("#!/usr/bin/env bash\n"
           "/home/ygeorgi/build/rotd/multi")
MCFLUX = ("#!/usr/bin/env bash\n"
          "/home/ygeorgi/build/rotd/mc_flux")
CONV_MULTI = ("#!/usr/bin/env bash\n"
              "/home/ygeorgi/build/rotd/mc_flux")
TST_CHECK = ("#!/usr/bin/env bash\n"
             "/home/ygeorgi/build/rotd/tst_check")

# Thermo
THERMP = ("#!/usr/bin/env bash\n"
          "THERMP FILE")
PAC99 = ("#!/usr/bin/env bash\n"
         "PACC << EOF\n"
         "FORMULA\n"
         "EOF")

# dsarrfit (likely not needed, as calling external not done)
DSARRFIT = ("#!/usr/bin/env bash\n"
            "dsarrfit.x_cfg")

# Electronic structure
G09 = ("#!/usr/bin/env bash\n"
       "g09 run.inp run.out >> stdout.log &> stderr.log")
PSI4 = ("#!/usr/bin/env bash\n"
        "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")
MOLPRO = ("#!/usr/bin/env bash\n"
          "molpro -n 4 run.inp -o run.out "
          "--nouse-logfile --no-xml-output >> "
          "stdout.log &> stderr.log")
MOLPRO_MPPX = ("#!/usr/bin/env bash\n"
               "molpro --mppx -n 4 run.inp -o run.out "
               "--nouse-logfile --no-xml-output >> "
               "stdout.log &> stderr.log")
