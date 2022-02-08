"""
Executes the automation part of 1DMin
"""

import ioformat
from mechlib.amech_io import output_path
from mechlib.amech_io.printer import obj, info_message
from mechroutines.trans._routines import run_onedmin
from mechroutines.trans._routines import build_transport_file


def run_tsk(tsk, spc_queue,
            spc_dct, thy_dct, pes_model_dct,
            etrans_keyword_dct,
            run_prefix, save_prefix, mdriver_path):
    """ Run the task
    """

    # Print the head of the task
    obj('vspace')
    obj('line_dash')
    info_message('Task:', tsk, newline=1)
    info_message('Options for transport task:', newline=1)

    for key, val in etrans_keyword_dct.items():
        info_message(f'{key} = {val}')

    if tsk == 'onedmin':
        obj('vspace')
        obj('line_dash')
        info_message(
            'Obtaining LJ-Params using OneDMin', newline=1)
        for spc_name, _ in spc_queue:
            run_onedmin(spc_name,
                        spc_dct, thy_dct, etrans_keyword_dct,
                        run_prefix, save_prefix)
    elif tsk == 'write_transport':

        # Read and Process data to determine transport properties
        obj('vspace')
        obj('line_dash')
        info_message(
            'TASK: Reading and Processing data into transport parameters',
            newline=1)
        ckin_str = build_transport_file(
            spc_queue,
            spc_dct, thy_dct, pes_model_dct,
            etrans_keyword_dct,
            run_prefix, save_prefix)

        # Write the file
        ckin_path = output_path('CKIN', prefix=mdriver_path)
        info_message(f'Writing CHEMKIN transport file at {ckin_path}')
        ioformat.pathtools.write_file(ckin_str, ckin_path, 'trans.ckin')
