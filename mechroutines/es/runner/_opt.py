"""
  Functions that are useful derived sequences of run_job and read_job
"""

import elstruct
from mechlib.amech_io import printer as ioprinter
from mechroutines.es.runner._run import execute_job


def multi_stage_optimization(script_str, run_fs,
                             geo, spc_info, thy_info,
                             frozen_coords_lst,
                             overwrite=False,
                             saddle=False,
                             retryfail=False,
                             **kwargs):
    """ Run a series of optimizations that utilize varying constraint
        conditions at each state.

        :param script_str: Shell submission script for electronic structure job
        :type script_str: str
        :param geo: input molecular geometry or Z-Matrix
        :type geo: automol.geom object
        :param spc_info:
        :type spc_info:
        :param thy_info:
        :type thy_info:
        :param frozen_coords_lst: Z-matrix coordinate names to freeze in opts
        :type frozen_coords_lst: tuple(str)
        :param overwrite: overwrite existing input file with new one and rerun
        :type overwrite: bool
        :param saddle: perform saddle-point optimization
        :type saddle: bool
        :param retryfail: re-run the job if failed job found in RUN filesys
        :type retryfail: bool
    """

    for idx, coords in enumerate(frozen_coords_lst):

        ioprinter.info_message(
            'Stage {} success beginning stage two'.format(idx+1))
        success, ret = execute_job(
            job=elstruct.Job.OPTIMIZATION,
            script_str=script_str,
            run_fs=run_fs,
            geo=geo,
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
            geo = elstruct.reader.opt_zmatrix(inf_obj.prog, out_str)
            print('Success. Moving to next stage...\n')
            if idx+1 != len(frozen_coords_lst):
                print('Success. Moving to next stage...\n')
            else:
                print('Succes. Finished sequence')
        else:
            print('Failure. Exiting multi stage optimization...\n')
            break

    return success, ret
