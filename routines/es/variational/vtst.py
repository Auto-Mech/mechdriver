""" Run and Read the scans from VTST calculations
"""

import automol
from routines.es import scan


def run_vtst_scan(ts_zma, ts_formula, ts_info, ts_dct, spc_dct,
                  high_mul, grid1, grid2, dist_name,
                  multi_level, num_act_orb, num_act_elc,
                  multi_info, ini_thy_info, thy_info,
                  run_prefix, save_prefix, scn_run_fs, scn_save_fs,
                  opt_script_str, overwrite, update_guess, **opt_kwargs):
    """ Run the scan for VTST calculations
    """
    gradient = False
    hessian = True

    moldr.scan.run_multiref_rscan(
        formula=ts_formula,
        high_mul=high_mul,
        zma=ts_zma,
        spc_info=ts_info,
        multi_level=multi_level,
        dist_name=dist_name,
        grid1=grid1,
        grid2=grid2,
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        script_str=opt_script_str,
        overwrite=overwrite,
        update_guess=update_guess,
        gradient=gradient,
        hessian=hessian,
        num_act_elc=num_act_elc,
        num_act_orb=num_act_orb,
        **opt_kwargs
    )

    moldr.scan.save_scan(
        scn_run_fs=scn_run_fs,
        scn_save_fs=scn_save_fs,
        coo_names=[dist_name],
        gradient=gradient,
        hessian=hessian,
    )

    locs = [[dist_name], [grid1[0]]]
    # calculate and save the infinite seperation energy
    print('ts zma locs')
    print(locs)
    ts_zma = scn_save_fs.leaf.file.zmatrix.read(locs)
    rcts = ts_dct['reacs']
    spc_1_info = [spc_dct[rcts[0]]['ich'],
                  spc_dct[rcts[0]]['chg'],
                  spc_dct[rcts[0]]['mul']]
    spc_2_info = [spc_dct[rcts[1]]['ich'],
                  spc_dct[rcts[1]]['chg'],
                  spc_dct[rcts[1]]['mul']]

    inf_sep_ene = moldr.scan.infinite_separation_energy(
        spc_1_info, spc_2_info, ts_info, high_mul, ts_zma, ini_thy_info, thy_info,
        multi_info, run_prefix, save_prefix, scn_run_fs, scn_save_fs, locs,
        num_act_elc=num_act_elc,
        num_act_orb=num_act_orb)

    inf_locs = [[dist_name], [1000.]]
    scn_save_fs.leaf.create(inf_locs)
    scn_save_fs.leaf.file.energy.write(inf_sep_ene, inf_locs)

    geo = automol.zmatrix.geometry(ts_zma)
    zma = ts_zma
    final_dist = grid1[0]

    return geo, zma, final_dist
