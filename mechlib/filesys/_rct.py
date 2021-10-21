""" rct cnf fs because don't knwo where else to put it and avoid
    circular imports
"""

from mechanalyzer.inf import thy as tinfo
from mechlib.filesys._build import build_fs
from mechlib.filesys.mincnf import min_energy_conformer_locators


def rcts_cnf_fs(
        rct_infos, thy_info, run_prefix, save_prefix,
        hbond_cutoffs=None):
    """ set reactant filesystem stuff
    """

    rct_cnf_fs = ()
    for rct_info in rct_infos:

        mod_thy_info = tinfo.modify_orb_label(
            thy_info, rct_info)

        # Build filesys for ini thy info
        cnf_run_fs, cnf_save_fs = build_fs(
            run_prefix, save_prefix, 'CONFORMER',
            spc_locs=rct_info,
            thy_locs=mod_thy_info[1:])

        loc_info = min_energy_conformer_locators(
            cnf_save_fs, mod_thy_info,
            hbond_cutoffs=hbond_cutoffs)
        min_cnf_locs, min_cnf_path = loc_info
        # Create run fs if that directory has been deleted to run the jobs
        if any(min_cnf_locs):
            if cnf_run_fs is not None:
                cnf_run_fs[-1].create(min_cnf_locs)
            rct_cnf_fs += ((cnf_run_fs, cnf_save_fs,
                            min_cnf_locs, min_cnf_path),)
        else:
            rct_cnf_fs += (None,)

    return rct_cnf_fs
