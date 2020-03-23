""" ts drivers
"""

import automol
import elstruct

# New libs
from routines.es import conformer
from lib.runner import driver
from lib.runner import par as runpar


def sadpt_reference_geometry(spcdct, thy_info, ini_thy_info,
                             thy_save_fs, ini_save_fs,
                             cnf_run_fs, cnf_save_fs,
                             run_fs, dist_info=(), overwrite=False):
    """ determine what to use as the reference geometry for all future runs
    If ini_thy_info refers to geometry dictionary then use that,
    otherwise values are from a hierarchy of:
    running level of theory, input level of theory, inchis.
    From the hierarchy an optimization is performed followed by a check for
    an imaginary frequency and then a conformer file system is set up.
    """

    thy_save_fs = geo_fs[3]
    ini_thy_save_fs = ini_fs[1]
    run_fs = geo_fs[-1]
    print('initializing geometry')
    geo = None
    geo_init = None
    # Check to see if geometry should be obtained from dictionary
    spc_info = [spcdct['ich'], spcdct['chg'], spcdct['mul']]
    if 'input_geom' in ini_thy_level:
        geom_obj = spcdct['geoobj']
        geo_init = geom_obj
        overwrite = True
        print('found initial geometry from geometry dictionary')
    else:
        # Check to see if geo already exists at running_theory
        if thy_save_fs[0].file.geometry.exists():
            thy_path = thy_save_fs[0].path()
            print(
                'getting reference geometry from {}'.format(thy_path))
            geo = thy_save_fs[0].file.geometry.read()
            zma = thy_save_fs[0].file.zmatrix.read()
            # print('geo:',automol.geom.string(geo))
        if not geo:
            if ini_thy_save_fs:
                if ini_thy_save_fs[0].file.geometry.exists():
                    thy_path = ini_thy_save_fs[0].path()
                    print(
                        'getting reference geometry from {}'.format(thy_path))
                    zma_init = ini_thy_save_fs[0].file.zmatrix.read()
                    geo_init = ini_thy_save_fs[0].file.geometry.read()
    if not geo and geo_init:
        _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(
            *thy_level[0:2])
        driver.run_job(
            job='optimization',
            script_str=opt_script_str,
            run_fs=run_fs,
            geom=zma_init,
            spc_info=spc_info,
            thy_level=thy_info,
            saddle=True,
            overwrite=overwrite,
            **opt_kwargs,
            )
        opt_ret = driver.read_job(
            job='optimization',
            run_fs=run_fs,
            )
        if opt_ret is not None:
            inf_obj, _, out_str = opt_ret
            prog = inf_obj.prog
            method = inf_obj.method
            ene = elstruct.reader.energy(prog, method, out_str)
            geo = elstruct.reader.opt_geometry(prog, out_str)
            zma = elstruct.reader.opt_zmatrix(prog, out_str)
        if geo:

            print(" - Saving...")
            print(" - Save path: {}".format(thy_save_fs[0].path()))

            thy_save_fs[0].file.energy.write(ene)
            thy_save_fs[0].file.geometry.write(geo)
            thy_save_fs[0].file.zmatrix.write(zma)

            conformer.single_conformer(
                spc_info, thy_info, 
                thy_save_fs, cnf_run_fs, cnf_save_fs,
                overwrite, saddle=saddle, dist_info=dist_info)
    return geo
