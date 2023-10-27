""" es_runners
"""

import numpy
import automol
import elstruct
import autofile
from phydat import phycon
from mechlib import filesys
from mechlib.amech_io.printer import reading, info_message
from mechlib.amech_io.printer import debug_message, warning_message
from mechlib.amech_io.printer import save_geo, save_energy
from mechroutines.es import runner as es_runner
from mechroutines.es._routines import _util as util


def tau_sampling(zma, ref_ene, spc_info,
                 mod_thy_info,
                 tau_run_fs, tau_save_fs,
                 script_str, overwrite,
                 db_style='directory',
                 nsamp_par=(False, 3, 3, 1, 50, 50),
                 tors_names=(),
                 zrxn=None, resave=False,
                 **kwargs):
    """ Sample over torsions optimizing all other coordinates
    """

    if resave:
        save_tau(
            tau_run_fs=tau_run_fs,
            tau_save_fs=tau_save_fs,
            mod_thy_info=mod_thy_info,
            db_style=db_style
        )

    run_tau(
        zma=zma,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        tau_run_fs=tau_run_fs,
        tau_save_fs=tau_save_fs,
        script_str=script_str,
        overwrite=overwrite,
        nsamp_par=nsamp_par,
        tors_names=tors_names,
        zrxn=zrxn,
        **kwargs,
    )

    save_tau(
        tau_run_fs=tau_run_fs,
        tau_save_fs=tau_save_fs,
        mod_thy_info=mod_thy_info,
        db_style=db_style
    )

    info_message(
        'Assessing the convergence of the Monte Carlo Partition Function...')
    assess_pf_convergence(tau_save_fs, ref_ene)


def run_tau(zma, spc_info, thy_info,
            tau_run_fs, tau_save_fs, script_str, overwrite,
            nsamp_par=(False, 3, 3, 1, 50, 50),
            tors_names=(),
            zrxn=None, **kwargs):
    """ run sampling algorithm to find tau dependent geometries
    """

    # Set the filesystem objects
    tau_save_fs[0].create()

    # Check if the structure is consistent with the filesystem
    _check_vma(zma, tau_save_fs)

    # Set the samples
    nsamp, tors_range_dct = util.calc_nsamp(
        tors_names, nsamp_par, zma, zrxn=zrxn)
    nsamp0 = nsamp
    nsampd = util.calc_nsampd(tau_save_fs, tau_run_fs, rid=None)

    num_to_samp = nsamp - nsampd

    info_message(
        ' - Number of samples that have been currently run:', nsampd)
    info_message(' - Number of samples requested:', nsamp)

    if num_to_samp > 0:
        info_message(
            f'Running {num_to_samp} samples...', newline=1)
    samp_idx = 1

    # Set the filesystem objects
    inf_obj = autofile.schema.info_objects.tau_trunk(0, tors_range_dct)

    while True:
        nsamp = nsamp0 - nsampd

        # Break the while loop if enough sampls completed
        if nsamp <= 0:
            info_message(
                'Reached requested number of samples. ',
                'Tau sampling complete.')
            break

        samp_zma, = automol.zmat.samples(zma, 1, tors_range_dct)
        tid = autofile.schema.generate_new_tau_id()
        locs = [tid]

        tau_run_fs[-1].create(locs)
        tau_run_prefix = tau_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(tau_run_prefix)

        info_message(f"\nRun {samp_idx}/{num_to_samp}")
        samp_idx += 1

        info_message(
            'Generating sample Z-Matrix that does not have',
            'high intramolecular repulsion...')
        if automol.zmat.has_low_relative_repulsion_energy(samp_zma, zma):
            debug_message('ZMA fine.')
            es_runner.run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geo=samp_zma,
                spc_info=spc_info,
                thy_info=thy_info,
                saddle=bool(zrxn is not None),
                overwrite=overwrite,
                frozen_coordinates=tors_range_dct.keys(),
                **kwargs
            )
        else:
            warning_message('repulsive ZMA:')
            inp_str = elstruct.writer.optimization(
                geo=samp_zma,
                charge=spc_info[1],
                mult=spc_info[2],
                method=thy_info[1],
                basis=thy_info[2],
                prog=thy_info[0],
                orb_type=thy_info[3],
                mol_options=['nosym'],
                frozen_coordinates=tors_range_dct.keys(),
            )
            tau_run_fs[-1].file.geometry_input.write(inp_str, locs)
            warning_message(
                'geometry for bad ZMA at', tau_run_fs[-1].path(locs))

        if tau_save_fs[0].file.info.exists():
            inf_obj_s = tau_save_fs[0].file.info.read()
            nsampd = inf_obj_s.nsamp
        elif tau_run_fs[0].file.info.exists():
            inf_obj_r = tau_run_fs[0].file.info.read()
            nsampd = inf_obj_r.nsamp
        nsampd += 1
        inf_obj.nsamp = nsampd
        tau_save_fs[0].file.info.write(inf_obj)
        tau_run_fs[0].file.info.write(inf_obj)


def save_tau(tau_run_fs, tau_save_fs, mod_thy_info, db_style='directory'):
    """ save the tau dependent geometries that have been found so far
    """

    if db_style == 'jsondb':
        saved_locs = tau_save_fs[-1].json_existing()
        saved_geos = tau_save_fs[-1].json.geometry.read_all(saved_locs)
    elif db_style == 'directory':
        saved_geos = [tau_save_fs[-1].file.geometry.read(locs)
                      for locs in tau_save_fs[-1].existing()]
    if not tau_run_fs[0].exists():
        info_message("No tau geometries to save. Skipping...")
    else:
        if db_style == 'jsondb':
            save_info = [[], [], [], [], []]
            sp_save_info = [[], [], [], [], []]
        for locs in tau_run_fs[-1].existing():
            run_path = tau_run_fs[-1].path(locs)
            run_fs = autofile.fs.run(run_path)
            save_path = tau_save_fs[-1].root.path()

            reading("tau run", run_path)

            success, ret = es_runner.read_job(
                job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if success:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)

                geo = elstruct.reader.opt_geometry(prog, out_str)
                if db_style == 'directory':
                    save_geo(save_path)
                    tau_save_fs[-1].create(locs)
                    tau_save_fs[-1].file.geometry_info.write(inf_obj, locs)
                    tau_save_fs[-1].file.geometry_input.write(inp_str, locs)
                    tau_save_fs[-1].file.energy.write(ene, locs)
                    tau_save_fs[-1].file.geometry.write(geo, locs)

                    # Saving the energy to a SP filesystem
                    save_path = tau_save_fs[-1].path(locs)
                    save_energy(save_path)
                    sp_save_fs = autofile.fs.single_point(save_path)
                    sp_save_fs[-1].create(mod_thy_info[1:4])
                    sp_save_fs[-1].file.input.write(inp_str, mod_thy_info[1:4])
                    sp_save_fs[-1].file.info.write(inf_obj, mod_thy_info[1:4])
                    sp_save_fs[-1].file.energy.write(ene, mod_thy_info[1:4])
                elif db_style == 'jsondb':
                    info_message(
                        " - Saving energy and geo of unique geometry...")
                    save_info[0].append(locs)
                    save_info[1].append(inf_obj)
                    save_info[2].append(inp_str)
                    save_info[3].append(ene)
                    save_info[4].append(geo)

                    sp_save_fs = autofile.fs.single_point(
                        save_path, json_layer=locs)
                    sp_save_info[0].append(sp_save_fs)
                    sp_save_info[1].append(mod_thy_info[1:4])
                    sp_save_info[2].append(inp_str)
                    sp_save_info[3].append(inf_obj)
                    sp_save_info[4].append(ene)
                    sp_save_fs[-1].json.input.write(inp_str, mod_thy_info[1:4])
                    sp_save_fs[-1].json.info.write(inf_obj, mod_thy_info[1:4])
                    sp_save_fs[-1].json.energy.write(ene, mod_thy_info[1:4])

                saved_geos.append(geo)

        print('\nWriting the geometries and energies into JSON file...')
        if db_style == 'jsondb':
            tau_save_fs[-1].json_create()
            tau_save_fs[-1].json.geometry_info.write_all(
                save_info[1], save_info[0])
            tau_save_fs[-1].json.geometry_input.write_all(
                save_info[2], save_info[0])
            tau_save_fs[-1].json.energy.write_all(
                save_info[3], save_info[0])
            tau_save_fs[-1].json.geometry.write_all(
                save_info[4], save_info[0])

            for i, sp_save_fs_i in enumerate(sp_save_info[0]):
                sp_save_fs_i[-1].json.input.write(
                    sp_save_info[2][i], sp_save_info[1][i])
                sp_save_fs_i[-1].json.info.write(
                    sp_save_info[3][i], sp_save_info[1][i])
                sp_save_fs_i[-1].json.energy.write(
                    sp_save_info[4][i], sp_save_info[1][i])

        # update the tau trajectory file
        filesys.mincnf.traj_sort(tau_save_fs, mod_thy_info)


def assess_pf_convergence(tau_save_fs, ref_ene,
                          temps=(300., 500., 750., 1000., 1500.)):
    """ Determine how much the partition function has converged
    """

    # Calculate sigma values at various temperatures for the PF
    for temp in temps:
        sumq = 0.
        sum2 = 0.
        idx = 0
        debug_message('integral convergence for T = ', temp)
        inf_obj_s = tau_save_fs[0].file.info.read()
        nsamp = inf_obj_s.nsamp
        saved_locs = tau_save_fs[-1].json_existing()
        ratio = len(saved_locs) / float(nsamp)
        for locs in saved_locs:
            idx += 1
            ene = tau_save_fs[-1].json.energy.read(locs)
            ene = (ene - ref_ene) * phycon.EH2KCAL
            tmp = numpy.exp(-ene*349.7/(0.695*temp))
            sumq = sumq + tmp
            sum2 = sum2 + tmp**2
            sigma = numpy.sqrt(
                (abs(sum2/float(idx)-(sumq/float(idx))**2))/float(idx))
            debug_message(
                sumq/float(idx), sigma, 100.*sigma*float(idx)/sumq, idx)
        info_message('Ratio of good to sampled geometries: ', ratio)


def _check_vma(zma, tau_save_fs):
    """ Assess of the vma matches the zma used to sample.
        Write the vma if needed.
    """

    vma = automol.zmat.vmatrix(zma)
    if tau_save_fs[0].file.vmatrix.exists():
        existing_vma = tau_save_fs[0].file.vmatrix.read()
        assert vma == existing_vma
    else:
        tau_save_fs[0].file.vmatrix.write(vma)
