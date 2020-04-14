""" Tests tasks associated with conformer tasks including
    sampling, energy, gradients, and hessians
"""

# SCRIPT_STR = (
#     "#!/bin/bash\n"
#     "python -u automech.py $CWD >& run.log & disown %1"
# )

SCRIPT_STR = (
    '#!/bin/bash\n'
    'export CWD=/lcrc/project/CMRP/runs/data/abst2\n'
    'export AMECHDIR=/lcrc/project/CMRP/amech\n'
    'export OUTPUTFILE=test.log\n'
    'ssh -n b447 " conda activate amech-env-mini                                                ;\n'
    '               . $AMECHDIR/fake-install.sh                                                 ;\n'
    '               cd $CWD                                                                     ;\n'
    '               module load gaussian/09-e.01                                                ;\n' 
    '               module load molpro/2015.1_170920                                            ;\n'
    '               python -u $AMECHDIR/moldriver/bin/automech.py $CWD >& $OUTPUTFILE &          \n'
    '               disown %1 "                                                                  \n'
) 


def test__ts():
    """ test all the ts searching  
    """


