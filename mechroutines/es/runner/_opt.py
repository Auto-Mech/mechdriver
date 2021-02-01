"""
  Functions that are useful derived sequences of run_job and read_job
"""

import elstruct
from mechroutines.es.runner._run import execute_job


def multi_stage_optimization(script_str, run_fs,
                             inp_geom, spc_info, thy_info,
                             frozen_coord_lst,
                             overwrite=False,
                             saddle=False,
                             retryfail=False,
                             **kwargs):
    """ Run some optimization where a set of coordinates are fixed
    """

    # Initialize loop geom to the input
    geom = inp_geom

    for idx, coords in enumerate(frozen_coord_lst):

        ioprinter.info_message('Stage one success beginning stage two')
        success, ret = execute_job(
            job=elstruct.Job.OPTIMIZATION,
            script_str=script_str,
            run_fs=run_fs,
            geom=geom,
            spc_info=spc_info,
            thy_info=thy_info,
            overwrite=overwrite,
            frozen_coordinates=coords,
            saddle=saddle,
            retryfail=retryfail,
            **kwargs
        )

        if success:
            inf_obj, _, out_str = ret
            geom = elstruct.reader.opt_zmatrix(inf_obj.prog, out_str)
            print('Success. Moving to next stage...\n')
            if idx+1 != len(frozen_coord_lst):
                print('Success. Moving to next stage...\n')
            else:
                print('Succes. Finished sequence')
        else:
            print('Failure. Exiting multi stage optimization...\n')
            break

    return success, ret
