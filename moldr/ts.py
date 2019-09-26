""" ts drivers
"""
from qcelemental import constants as qcc
import automol
import elstruct
import autofile
import moldr
import scripts.es

WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

def ts_conformer_sampling(
        spc_info, geo, zma, tors_names, thy_level, thy_save_fs, cnf_run_fs, cnf_save_fs,
        script_str, overwrite, saddle=True, nsamp_par=(False, 3, 3, 1, 50, 50),
        **kwargs):
    """ Find the minimum energy conformer by optimizing from nsamp random
    initial torsional states
    """
    tors_ranges = automol.zmatrix.torsional_sampling_ranges(
        zma, tors_names)
    tors_range_dct = dict(zip(tors_names, tors_ranges))
    ntaudof = len(tors_names)
    nsamp = moldr.util.nsamp_init(nsamp_par, ntaudof)

    ts_save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
    )

    moldr.conformer.run_conformers(
        zma=zma,
        spc_info=spc_info,
        thy_level=thy_level,
        nsamp=nsamp,
        tors_range_dct=tors_range_dct,
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        script_str=script_str,
        overwrite=overwrite,
        **kwargs,
    )

    ts_save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
    )

    # save information about the minimum energy conformer in top directory
    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        zma = cnf_save_fs.leaf.file.zmatrix.read(min_cnf_locs)

        thy_save_fs.leaf.file.geometry.write(geo, thy_level[1:4])
        thy_save_fs.leaf.file.zmatrix.write(zma, thy_level[1:4])


def ts_save_conformers(cnf_run_fs, cnf_save_fs):
    """ save the conformers that have been found so far
    """

    locs_lst = cnf_save_fs.leaf.existing()
    seen_geos = [cnf_save_fs.leaf.file.geometry.read(locs)
                 for locs in locs_lst]
    seen_enes = [cnf_save_fs.leaf.file.energy.read(locs)
                 for locs in locs_lst]

    if not cnf_run_fs.trunk.exists():
        print("No conformers to save. Skipping...")
    else:
        for locs in cnf_run_fs.leaf.existing():
            run_path = cnf_run_fs.leaf.path(locs)
            run_fs = autofile.fs.run(run_path)
            print("Reading from conformer run at {}".format(run_path))

            ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                zma = elstruct.reader.opt_zmatrix(prog, out_str)

                unique = True

                for idx, geoi in enumerate(seen_geos):
                    enei = seen_enes[idx]
                    etol = 1.e-6
                    if automol.geom.almost_equal_coulomb_spectrum(
                            geo, geoi, rtol=1e-2):
                        if abs(ene-enei) < etol:
                            unique = False
                if not unique:
                    print(" - Geometry is not unique. Skipping...")
                else:
                    save_path = cnf_save_fs.leaf.path(locs)
                    print(" - Geometry is unique. Saving...")
                    print(" - Save path: {}".format(save_path))

                    cnf_save_fs.leaf.create(locs)
                    cnf_save_fs.leaf.file.geometry_info.write(
                        inf_obj, locs)
                    cnf_save_fs.leaf.file.geometry_input.write(
                        inp_str, locs)
                    cnf_save_fs.leaf.file.energy.write(ene, locs)
                    cnf_save_fs.leaf.file.geometry.write(geo, locs)
                    cnf_save_fs.leaf.file.zmatrix.write(zma, locs)

                seen_geos.append(geo)
                seen_enes.append(ene)

        # update the conformer trajectory file
        locs_lst = cnf_save_fs.leaf.existing()
        if locs_lst:
            enes = [cnf_save_fs.leaf.file.energy.read(locs)
                    for locs in locs_lst]
            geos = [cnf_save_fs.leaf.file.geometry.read(locs)
                    for locs in locs_lst]

            traj = []
            for ene, geo in sorted(zip(enes, geos), key=lambda x: x[0]):
                comment = 'energy: {:>15.10f}'.format(ene)
                traj.append((comment, geo))

            traj_path = cnf_save_fs.trunk.file.trajectory.path()
            print("Updating conformer trajectory file at {}".format(traj_path))
            cnf_save_fs.trunk.file.trajectory.write(traj)


def reference_geometry(
    spcdct, thy_level, ini_thy_level, fs, ini_fs, dist_info=[], overwrite=False):
    """ determine what to use as the reference geometry for all future runs
    If ini_thy_info refers to geometry dictionary then use that,
    otherwise values are from a hierarchy of:
    running level of theory, input level of theory, inchis.
    From the hierarchy an optimization is performed followed by a check for
    an imaginary frequency and then a conformer file system is set up.
    """

    thy_save_fs = fs[3]
    ini_thy_save_fs = ini_fs[1]
    run_fs = fs[-1]
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
        if thy_save_fs.trunk.file.geometry.exists():
            thy_path = thy_save_fs.trunk.path()
            print('getting reference geometry from {}'.format(thy_path))
            geo = thy_save_fs.trunk.file.geometry.read()
            zma = thy_save_fs.trunk.file.zmatrix.read()
        if not geo:
            if ini_thy_save_fs:
                if ini_thy_save_fs.trunk.file.geometry.exists():
                    thy_path = ini_thy_save_fs.trunk.path()
                    print('getting reference geometry from {}'.format(thy_path))
                    zma_init = ini_thy_save_fs.trunk.file.zmatrix.read()
                    geo_init = ini_thy_save_fs.trunk.file.geometry.read()
    if not geo and geo_init:
        _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(*thy_level[0:2])
        moldr.driver.run_job(
            job='optimization',
            script_str=opt_script_str,
            run_fs=run_fs,
            geom=zma_init,
            spc_info=spc_info,
            thy_level=thy_level,
            saddle=True,
            overwrite=overwrite,
            **opt_kwargs,
            )
        opt_ret = moldr.driver.read_job(
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
            print(" - Save path: {}".format(thy_save_fs.trunk.path()))

            thy_save_fs.trunk.file.energy.write(ene)
            thy_save_fs.trunk.file.geometry.write(geo)
            thy_save_fs.trunk.file.zmatrix.write(zma)

            scripts.es.run_single_conformer(
                spc_info, thy_level, fs,
                overwrite, True, dist_info)
    return geo
