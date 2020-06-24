"""
  Submit multiple automech jobs at once
"""

import os

CWD = os.getcwd()

RUN_DCT = mrun.submit()

for run, vals in RUN_DCT.items():
    # make the run dir
    # copy input files over
    # add the pes,chnls to run.dat
    # write the submission file
    # submit job if wanted

