""" rct cnf fs because don't knwo where else to put it and avoid
    circular imports
"""

from mechanalyzer.inf import thy as tinfo
from mechlib.filesys._build import build_fs
from mechlib.filesys.mincnf import min_energy_conformer_locators


def rcts_cnf_fs(rct_infos, thy_dct, es_keyword_dct, run_prefix, save_prefix):
    """ set reactant filesystem stuff
    """

    ini_method_dct = thy_dct.get(es_keyword_dct['inplvl'])
    ini_thy_info = tinfo.from_dct(ini_method_dct)

    rct_cnf_fs = ()
    for rct_info in rct_infos:

        mod_ini_thy_info = tinfo.modify_orb_label(
            ini_thy_info, rct_info)

        # Build filesys for ini thy info
        ini_cnf_run_fs, ini_cnf_save_fs = build_fs(
            run_prefix, save_prefix, 'CONFORMER',
            spc_locs=rct_info,
            thy_locs=mod_ini_thy_info[1:])

        ini_loc_info = min_energy_conformer_locators(
            ini_cnf_save_fs, mod_ini_thy_info)
        ini_min_cnf_locs, ini_min_cnf_path = ini_loc_info

        # Create run fs if that directory has been deleted to run the jobs
        if any(ini_min_cnf_locs):
            ini_cnf_run_fs[-1].create(ini_min_cnf_locs)
            rct_cnf_fs += ((ini_cnf_run_fs, ini_cnf_save_fs,
                            ini_min_cnf_locs, ini_min_cnf_path),)
        else:
            rct_cnf_fs += ((),)

    return rct_cnf_fs
