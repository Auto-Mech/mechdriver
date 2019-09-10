""" drivers for conformer sampling 
"""
from qcelemental import constants as qcc
import automol
import elstruct
import autofile
import moldr

WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

def conformer_sampling(
        spc_info, thy_level, thy_save_fs, cnf_run_fs, cnf_save_fs,
        script_str, overwrite, saddle=False, nsamp_par=(False, 3, 3, 1, 50, 50), tors_names='',
        **kwargs):
    """ Find the minimum energy conformer by optimizing from nsamp random
    initial torsional states
    """

    geo = thy_save_fs.leaf.file.geometry.read(thy_level[1:4])
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    zma = automol.geom.zmatrix(geo)
    tors_ranges = automol.zmatrix.torsional_sampling_ranges(
        zma, tors_names)
    tors_range_dct = dict(zip(tors_names, tors_ranges))
    ich = spc_info[0]
    gra = automol.inchi.graph(ich)
    ntaudof = len(automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
    nsamp = moldr.util.nsamp_init(nsamp_par, ntaudof)

    save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
    )
#    print('zma test in conformer sampling')
#    print(automol.zmatrix.string(zma))

    run_conformers(
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

    save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
    )

    # save information about the minimum energy conformer in top directory
    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        zma = cnf_save_fs.leaf.file.zmatrix.read(min_cnf_locs)
        assert automol.zmatrix.almost_equal(zma, automol.geom.zmatrix(geo))
        thy_save_fs.leaf.file.geometry.write(geo, thy_level[1:4])
        thy_save_fs.leaf.file.zmatrix.write(zma, thy_level[1:4])


def run_conformers(
        zma, spc_info, thy_level, nsamp, tors_range_dct,
        cnf_run_fs, cnf_save_fs, script_str, overwrite,
        **kwargs):
    """ run sampling algorithm to find conformers
    """
    if not tors_range_dct:
        print("No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    print('Number of samples requested:', nsamp)

    cnf_save_fs.trunk.create()

    vma = automol.zmatrix.var_(zma)
    if cnf_save_fs.trunk.file.vmatrix.exists():
        existing_vma = cnf_save_fs.trunk.file.vmatrix.read()
        assert vma == existing_vma
    cnf_save_fs.trunk.file.vmatrix.write(vma)
    idx = 0
    nsamp0 = nsamp
    inf_obj = autofile.system.info.conformer_trunk(0, tors_range_dct)
    while True:
        if cnf_save_fs.trunk.file.info.exists():
            inf_obj_s = cnf_save_fs.trunk.file.info.read()
            nsampd = inf_obj_s.nsamp
        elif cnf_run_fs.trunk.file.info.exists():
            inf_obj_r = cnf_run_fs.trunk.file.info.read()
            nsampd = inf_obj_r.nsamp
        else:
            nsampd = 0

        nsamp = nsamp0 - nsampd
        if nsamp <= 0:
            print('Reached requested number of samples. '
                  'Conformer search complete.')
            break
        else:
            print("    New nsamp requested is {:d}.".format(nsamp))

            samp_zma, = automol.zmatrix.samples(zma, 1, tors_range_dct)
            cid = autofile.system.generate_new_conformer_id()
            locs = [cid]

            cnf_run_fs.leaf.create(locs)
            cnf_run_path = cnf_run_fs.leaf.path(locs)
            run_fs = autofile.fs.run(cnf_run_path)
#            print('run_conformers test:', cid, cnf_run_path)

            idx += 1
            print("Run {}/{}".format(idx, nsamp0))
            moldr.driver.run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geom=samp_zma,
                spc_info=spc_info,
                thy_level=thy_level,
                overwrite=overwrite,
                **kwargs
            )

            nsampd += 1
            inf_obj.nsamp = nsampd
            cnf_save_fs.trunk.file.info.write(inf_obj)
            cnf_run_fs.trunk.file.info.write(inf_obj)


def save_conformers(cnf_run_fs, cnf_save_fs):
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
            cnf_run_path = cnf_run_fs.leaf.path(locs)
            run_fs = autofile.fs.run(cnf_run_path)
            print("Reading from conformer run at {}".format(cnf_run_path))

            ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                zma = automol.geom.zmatrix(geo)
                gra = automol.geom.graph(geo)

                if len(automol.graph.connected_components(gra)) > 1:
                    print(" - Geometry is disconnected.. Skipping...")
                else:
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


def run_conformer_gradients(
        spc_info, thy_level, cnf_run_fs, cnf_save_fs,
        script_str, overwrite, **kwargs):
    """ Determine the gradient for each of the conformers
    """

    cnf_locs_lst = cnf_save_fs.leaf.existing()
    for locs in cnf_locs_lst:
        cnf_run_path = cnf_run_fs.leaf.path(locs)
        cnf_save_path = cnf_save_fs.leaf.path(locs)
        geo = cnf_save_fs.leaf.file.geometry.read(locs)
        run_fs = autofile.fs.run(cnf_run_path)

        print('Running conformer gradient')
        moldr.driver.run_job(
            job='gradient',
            script_str=script_str,
            run_fs=cnf_run_fs,
            geom=geo,
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            **kwargs,
        )

        ret = moldr.driver.read_job(
            job='gradient',
            run_fs=run_fs,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            print(" - Reading gradient from output...")
            grad = elstruct.reader.gradient(inf_obj.prog, out_str)

            print(" - Saving gradient...")
            print(" - Save path: {}".format(cnf_save_path))
            cnf_save_fs.leaf.file.gradient_info.write(inf_obj, locs)
            cnf_save_fs.leaf.file.gradient_input.write(inp_str, locs)
            cnf_save_fs.leaf.file.gradient.write(grad, locs)


def run_conformer_hessians(
        spc_info, thy_level, cnf_run_fs, cnf_save_fs,
        script_str, overwrite, **kwargs):
    """ Determine the hessian for each of the conformers
    """

    cnf_locs_lst = cnf_save_fs.leaf.existing()
    for locs in cnf_locs_lst:
        cnf_run_path = cnf_run_fs.leaf.path(locs)
        cnf_save_path = cnf_save_fs.leaf.path(locs)
        geo = cnf_save_fs.leaf.file.geometry.read(locs)
        run_fs = autofile.fs.run(cnf_run_path)

        print('Running conformer hessian')
        moldr.driver.run_job(
            job='hessian',
            script_str=script_str,
            run_fs=cnf_run_fs,
            geom=geo,
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            **kwargs,
        )

        ret = moldr.driver.read_job(
            job='hessian',
            run_fs=run_fs,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            if automol.geom.is_atom(geo):
                freqs = ()
            else:
                print(" - Reading hessian from output...")
                hess = elstruct.reader.hessian(inf_obj.prog, out_str)
                freqs = elstruct.util.harmonic_frequencies(
                    geo, hess, project=False)
#                print('Conformer Freqs test')
#                print(freqs)

                print(" - Saving hessian...")
                print(" - Save path: {}".format(cnf_save_path))
                cnf_save_fs.leaf.file.hessian_info.write(inf_obj, locs)
                cnf_save_fs.leaf.file.hessian_input.write(inp_str, locs)
                cnf_save_fs.leaf.file.hessian.write(hess, locs)
                cnf_save_fs.leaf.file.harmonic_frequencies.write(freqs, locs)

