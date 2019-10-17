""" test
"""

import read_dat

# Read the run parameters from a datafile
PARAMS = read_dat.params('params.dat')

for attr in dir(PARAMS):
    if '__' not in attr:
        print(attr, getattr(PARAMS, attr))
