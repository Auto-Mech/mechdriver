"""
Executes the automation part of 1DMin
"""

from routines.trans._routines import lj
from routines.trans._routines import build


def run_tsk(tsk, spc_queue,
            spc_dct,
            thy_dct, etrans_keyword_dct,
            run_prefix, save_prefix):
    """ Run the task
    """

    # Print the head of the task
    print(('\n------------------------------------------------' +
           '--------------------------------------'))
    print('\nTask:', tsk)
    print('\nOptions for transport task:')
    for key, val in etrans_keyword_dct.items():
        print('{} = {}'.format(key, val))
    print('')

    if tsk == 'onedmin':
        print(('\n\n------------------------------------------------' +
               '--------------------------------------'))
        print('\nObtaining LJ-Params using OneDMin')
        for spc_name, _ in spc_queue:
            lj.onedmin(spc_name,
                       spc_dct, thy_dct, etrans_keyword_dct,
                       run_prefix, save_prefix)
    elif tsk == 'write_transport':
        print(('\n\n------------------------------------------------' +
               '--------------------------------------'))
        print('\nWriting the CHEMKIN transport file')
        build.collate_properties(spc_queue,
                                 spc_dct, thy_dct, etrans_keyword_dct,
                                 save_prefix)
