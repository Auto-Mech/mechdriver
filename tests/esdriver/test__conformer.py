""" Tests tasks associated with conformer tasks including
    sampling, energy, gradients, and hessians
"""

SCRIPT_STR = (
    "#!/bin/bash\n"
    "python -u automech.py $CWD >& run.log & disown %1"
)


