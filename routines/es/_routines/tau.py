""" es_runners
"""

# import itertools
import numpy
import automol
import elstruct
import autofile
from routines.es._routines import _util as util
from routines.es import runner as es_runner
from lib import filesys
from lib.phydat import phycon


def tau_sampling(zma, ref_ene, spc_info, tors_name_grps, nsamp_par,
                 mod_thy_info,
                 tau_run_fs, tau_save_fs,
                 script_str, overwrite,
                 saddle=False, **opt_kwargs):
    """ Sample over torsions optimizing all other coordinates
    """

    # Read the geometry from the initial filesystem and set sampling
    tors_ranges = automol.zmatrix.torsional_sampling_ranges(tors_name_grps)
    tors_range_dct = dict(zip(
        tuple(grp[0] for grp in tors_name_grps), tors_ranges))
    gra = automol.zmatrix.graph(zma)
    ntaudof = len(automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
    nsamp = util.nsamp_init(nsamp_par, ntaudof)
    # Run through tau sampling process
    save_tau(
        tau_run_fs=tau_run_fs,
        tau_save_fs=tau_save_fs,
        mod_thy_info=mod_thy_info
    )

    run_tau(
        zma=zma,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        nsamp=nsamp,
        tors_range_dct=tors_range_dct,
        tau_run_fs=tau_run_fs,
        tau_save_fs=tau_save_fs,
        script_str=script_str,
        overwrite=overwrite,
        saddle=saddle,
        **opt_kwargs,
    )

    save_tau(
        tau_run_fs=tau_run_fs,
        tau_save_fs=tau_save_fs,
        mod_thy_info=mod_thy_info
    )

    print('Assessing the convergence of the Monte Carlo Partition Function...')
    assess_pf_convergence(tau_save_fs, ref_ene)


def run_tau(zma, spc_info, thy_info, nsamp, tors_range_dct,
            tau_run_fs, tau_save_fs, script_str, overwrite,
            saddle, **kwargs):
    """ run sampling algorithm to find tau dependent geometries
    """
    if not tors_range_dct:
        print("No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    tau_save_fs[0].create()

    vma = automol.zmatrix.var_(zma)
    if tau_save_fs[0].file.vmatrix.exists():
        existing_vma = tau_save_fs[0].file.vmatrix.read()
        assert vma == existing_vma
    tau_save_fs[0].file.vmatrix.write(vma)
    idx = 0
    nsamp0 = nsamp
    inf_obj = autofile.schema.info_objects.tau_trunk(0, tors_range_dct)
    while True:
        if tau_save_fs[0].file.info.exists():
            inf_obj_s = tau_save_fs[0].file.info.read()
            nsampd = inf_obj_s.nsamp
        elif tau_run_fs[0].file.info.exists():
            inf_obj_r = tau_run_fs[0].file.info.read()
            nsampd = inf_obj_r.nsamp
        else:
            nsampd = 0

        nsamp = nsamp0 - nsampd
        if nsamp <= 0:
            print('Reached requested number of samples. '
                  'Tau sampling complete.')
            break

        print("    New nsamp is {:d}.".format(nsamp))

        samp_zma, = automol.zmatrix.samples(zma, 1, tors_range_dct)
        tid = autofile.schema.generate_new_tau_id()
        locs = [tid]

        tau_run_fs[-1].create(locs)
        tau_run_prefix = tau_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(tau_run_prefix)

        idx += 1
        print("Run {}/{}".format(idx, nsamp0))

        print('\nChecking if ZMA has high repulsion...')
        if automol.intmol.low_repulsion_struct(zma, samp_zma):
            print('ZMA fine.')
            es_runner.run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geom=samp_zma,
                spc_info=spc_info,
                thy_info=thy_info,
                saddle=saddle,
                overwrite=overwrite,
                frozen_coordinates=tors_range_dct.keys(),
                **kwargs
            )
        else:
            print('repulsive ZMA:')
            inp_str = elstruct.writer.optimization(
                geom=samp_zma,
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
            print('geometry for bad ZMA at', tau_run_fs[-1].path(locs))

        # nsampd += 1
        # inf_obj.nsamp = nsampd
        # tau_save_fs[0].file.info.write(inf_obj)
        # tau_run_fs[0].file.info.write(inf_obj)

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


def save_tau(tau_run_fs, tau_save_fs, mod_thy_info):
    """ save the tau dependent geometries that have been found so far
    """
    #saved_geos = [tau_save_fs[-1].json.geometry.read(locs)
     #             for locs in tau_save_fs[-1].json_existing()]
    saved_locs = tau_save_fs[-1].existing() 
    saved_geos = tau_save_fs[-1].json.geometry.read_all(saved_locs)
    if not tau_run_fs[0].exists():
        print("No tau geometries to save. Skipping...")
    else:
        save_info = [[], [], [], [], []]
        sp_save_info = [[], [], [], [], []]
        for locs in tau_run_fs[-1].existing():
            run_path = tau_run_fs[-1].path(locs)
            run_fs = autofile.fs.run(run_path)
            save_path = tau_save_fs[-1].root.path()

            print("Reading from tau run at {}".format(run_path))

            success, ret = es_runner.read_job(
                job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if success:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)

                geo = elstruct.reader.opt_geometry(prog, out_str)
#
#                save_path = tau_save_fs[-1].path(locs)
#                print(" - Saving...")
#                print(" - Save path: {}".format(save_path))
                print(" - Saving...")
                print(" - Save path: {}".format(save_path))

#                tau_save_fs[-1].create(locs)
#                tau_save_fs[-1].file.geometry_info.write(inf_obj, locs)
#                tau_save_fs[-1].file.geometry_input.write(inp_str, locs)
#                tau_save_fs[-1].file.energy.write(ene, locs)
#                tau_save_fs[-1].file.geometry.write(geo, locs)
#                
                save_info[0].append(locs)        
                save_info[1].append(inf_obj)        
                save_info[2].append(inp_str)        
                save_info[3].append(ene)        
                save_info[4].append(geo)        
              #  tau_save_fs[-1].json.geometry_info.write(inf_obj, locs)
              #  tau_save_fs[-1].json.geometry_input.write(inp_str, locs)
              #  tau_save_fs[-1].json.energy.write(ene, locs)
              #  tau_save_fs[-1].json.geometry.write(geo, locs)

                # Saving the energy to a SP filesystem

                # Saving the energy to a SP filesystem
                print(" - Saving energy of unique geometry...")
#                sp_save_fs = autofile.fs.single_point(save_path)
#                sp_save_fs[-1].create(mod_thy_info[1:4])
#                sp_save_fs[-1].file.input.write(inp_str, mod_thy_info[1:4])
#                sp_save_fs[-1].file.info.write(inf_obj, mod_thy_info[1:4])
#                sp_save_fs[-1].file.energy.write(ene, mod_thy_info[1:4])
                sp_save_fs = autofile.fs.single_point(save_path, json_layer=locs)
                sp_save_info[0].append(sp_save_fs)
                sp_save_info[1].append(mod_thy_info[1:4])
                sp_save_info[2].append(inp_str)
                sp_save_info[3].append(inf_obj)
                sp_save_info[4].append(ene)
                sp_save_fs[-1].json.input.write(inp_str, mod_thy_info[1:4])
                sp_save_fs[-1].json.info.write(inf_obj, mod_thy_info[1:4])
                sp_save_fs[-1].json.energy.write(ene, mod_thy_info[1:4])

                saved_geos.append(geo)

        tau_save_fs[-1].json_create()
        tau_save_fs[-1].json.geometry_info.write_all(save_info[1], save_info[0])
        tau_save_fs[-1].json.geometry_input.write_all(save_info[2], save_info[0])
        tau_save_fs[-1].json.energy.write_all(save_info[3], save_info[0])
        tau_save_fs[-1].json.geometry.write_all(save_info[4], save_info[0])

        for i, sp_save_fs_i in enumerate(sp_save_info[0]):
            sp_save_fs_i[-1].json.input.write(sp_save_info[2][i], sp_save_info[1][i])
            sp_save_fs_i[-1].json.info.write(sp_save_info[3][i], sp_save_info[1][i])
            sp_save_fs_i[-1].json.energy.write(sp_save_info[4][i], sp_save_info[1][i])

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
        print('integral convergence for T = ', temp)
        for locs in tau_save_fs[-1].existing():
            idx += 1
            ene = tau_save_fs[-1].file.energy.read(locs)
            ene = (ene - ref_ene) * phycon.EH2KCAL
            tmp = numpy.exp(-ene*349.7/(0.695*temp))
            sumq = sumq + tmp
            sum2 = sum2 + tmp**2
            sigma = numpy.sqrt(
                (abs(sum2/float(idx)-(sumq/float(idx))**2))/float(idx))
            print(sumq/float(idx), sigma, 100.*sigma*float(idx)/sumq, idx)
        inf_obj_s = tau_save_fs[0].file.info.read()
        nsamp = inf_obj_s.nsamp
        saved_locs = tau_save_fs[-1].existing()
        ratio = len(saved_locs) / float(nsamp)
        print('ratio of good to sampled geometries', ratio)
