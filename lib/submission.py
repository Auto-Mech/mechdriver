""" Run bash scripts
"""

import os
# import sys
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


# Build a dictionary of submission scripts
# SUB_DIR = 'inp/sub_scripts'
# def read_inp_scripts(job_path):
#     """ Read in the scripts the user submits as a strings
#     """
#
#     # Initialize dct with that from standard defined below
#     sub_dct = DEFAULT_SCRIPT_DCT
#
#     # Obtain a list of any submission script files from input dir
#     sub_path = os.path.join(job_path, SUB_DIR)
#     if os.path.exists(sub_path):
#         sub_files = os.listdir(sub_path)
#
#     # Build dct if any sub files found
#     if sub_files:
#         print('Found directory with submission scripts.')
#         for sub_file in sub_files:
#             name = sfile.replace('.sh', '')
#             if name in sub_dct:
#                 sub_file_path = os.path.join(sub_path, sub_file)
#                 print('  Reading submission script for {} at path {}'.format(
#                     name, sub_file_path))
#                 with open(sub_file_path, 'r') as sfile:
#                     sub_dct[name] = script_file.read()
#             else:
#                 print('  Submission script unviable at path {}'.format(
#                     sub_file_path)
#                 print('  Allowed file names: ')
#                 for key in DEFAULT_SCRIPT_DCT:
#                     print('   ', key)
#                 sys.exit()
#     else:
#         print('No submmission script directory provided by the user')
#
#     return script_dct


# DCT OF DEFAULT SUBMISSION STRINGS
PROJROT = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "RPHt.exe >& /dev/null"
)
MESSPF = (
    "#!/usr/bin/env bash\n"
    "export OMP_NUM_THREADS=10\n"
    "ulimit -c 0\n"
    "messpf pf.inp pf.out >> stdout.log &> stderr.log"
)
MESSRATE = (
    "#!/usr/bin/env bash\n"
    "export OMP_NUM_THREADS=10\n"
    "ulimit -c 0\n"
    "mess mess.inp rate.out >> stdout.log &> stderr.log"
)
VARECOF = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/home/ygeorgi/build/rotd/multi >& varecof.out"
)
MCFLUX = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/home/ygeorgi/build/rotd/mc_flux"
)
TSTCHECK = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/home/ygeorgi/build/rotd/tst_check"
)
THERMP = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "THERMP FILE"
)
PAC99 = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "PACC << EOF\n"
    "FORMULA\n"
    "EOF"
)
DSARRFIT = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "dsarrfit.x_cfg"
)
G09 = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "g09 run.inp run.out >> stdout.log &> stderr.log"
)
PSI4 = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log"
)
MOLPRO = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "molpro -n 4 run.inp -o run.out "
    "--nouse-logfile --no-xml-output >> "
    "stdout.log &> stderr.log"
)
MOLPRO_MPPX = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "molpro --mppx -n 4 run.inp -o run.out "
    "--nouse-logfile --no-xml-output >> "
    "stdout.log &> stderr.log"
)
MOLPRO_MREF = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "molpro -n 8 run.inp -o run.out "
    "--nouse-logfile --no-xml-output >> "
    "stdout.log &> stderr.log"
)
MOLPRO_MREF_MPPX = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "molpro --mppx -n 12 run.inp -o run.out "
    "--nouse-logfile --no-xml-output >> "
    "stdout.log &> stderr.log"
)

DEFAULT_SCRIPT_DCT = {
    'projrot': PROJROT,
    'messpf': MESSPF,
    'messrate': MESSRATE,
    'varecof': VARECOF,
    'mcflux': MCFLUX,
    'tstchk': TSTCHECK,
    'thermp': THERMP,
    'pac99': PAC99,
    'dsarrfit': DSARRFIT,
    'gaussian09': G09,
    'molpro2015': MOLPRO,
    'molpro2015_mppx': MOLPRO_MPPX,
    'molpro2015_mr': MOLPRO_MREF,
    'molpro2015_mr_mppx': MOLPRO_MREF_MPPX
}
