"""
A simple script for running amech across multiple PESs
"""

import os
import sys
import subprocess
import argparse


RUN_TEMP = 'inp/run.temp'
EXE_TEMP = 'inp/exe.temp'
NEW_RUN = 'inp/run.dat'
NEW_EXE = 'exe.sh'

# Check for the existence of template files and run them
if os.path.exists(RUN_TEMP):
    with open(RUN_TEMP, 'r') as runfile:
        temp_run_str = runfile.read()
else:
    print('NO run.temp file in inp directory')
    sys.exit()
if os.path.exists(EXE_TEMP):
    with open(EXE_TEMP, 'r') as exefile:
        temp_exe_str = exefile.read()
else:
    print('NO exe.temp file in inp directory')
    sys.exit()

# Set all of the running stuff
cmd_line_parser = argparse.ArgumentParser()
cmd_line_parser.add_argument("pes_idx", help="Idx for pes to run")
cmd_line_parser.add_argument("chn_idxs", help="Idxs for channels to run (int: 1, or hyphen-delim list: 1-6)")
cmd_line_parser.add_argument("pes_model", help="Name of the pes model to use for the channel")
cmd_line_parser.add_argument("spc_model", help="Name of the spc model to use for the channel")
cmd_line_parser.add_argument("node", help="Name of the node to submit to")
cmd_line_parser.add_argument("log_file", help="Name of the log file for this run")
options = vars(cmd_line_parser.parse_args())

# Format channel list for run.dat
if '-' in options['chn_idxs']:
    [cstart, cend] = options['chn_idxs'].split('-')
    chn_lst = [str(num) for num in range(int(cstart), int(cend)+1)]
else:
    chn_lst = [options['chn_idxs']]

chn_str = ''
for i, num in enumerate(chn_lst):
    chn_str += '    {} {} ; {} {}'.format(options['pes_idx'], options['pes_model'], num, options['spc_model'])
    if i+1 != len(chn_lst):
        chn_str += '\n'

# Replace template strings 
new_run_str = temp_run_str.replace('RUN_LST', chn_str)
new_exe_str = temp_exe_str.replace('NODE', options['node']).replace('LOGFILE', options['log_file'])

# Write the new run.dat and exe.sh files
with open(NEW_RUN, 'w') as runfile:
    runfile.write(new_run_str)
with open(NEW_EXE, 'w') as exefile:
    exefile.write(new_exe_str)

# Execute the exe script
subprocess.check_call(['chmod', '+x', NEW_EXE])
subfile = open('submit.log', 'w')
subprocess.call(["./"+NEW_EXE], stdout=subfile, stderr=subprocess.STDOUT)
subfile.close()
print('Job submitted to node: '+options['node']+'\n')


