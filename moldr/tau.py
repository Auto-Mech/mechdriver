""" drivers
"""
import os
import numpy
from qcelemental import constants as qcc
import automol
import elstruct
import autofile
import moldr

WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')


def tau_sampling(
        spc_info, thy_level, thy_save_fs, tau_run_fs, tau_save_fs,
        script_str, overwrite, nsamp_par, **opt_kwargs):
    """ Sample over torsions optimizing all other coordinates
    """
#    thy_run_fs.leaf.create(thy_level)
#    thy_run_path = thy_run_fs.leaf.path(thy_level)

#    thy_save_fs.leaf.create(thy_level)
#    thy_save_path = thy_save_fs.leaf.path(thy_level)
#    tau_run_fs = autofile.fs.tau(run_prefix)
#    tau_save_fs = autofile.fs.tau(save_prefix)

    save_tau(
        tau_run_fs=tau_run_fs,
        tau_save_fs=tau_save_fs,
    )

    geo = thy_save_fs.leaf.file.geometry.read(thy_level[1:4])
    zma = automol.geom.zmatrix(geo)
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    tors_ranges = automol.zmatrix.torsional_sampling_ranges(
        zma, tors_names)
    tors_range_dct = dict(zip(tors_names, tors_ranges))
    ich = spc_info[0]
    gra = automol.inchi.graph(ich)
    ntaudof = len(automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
    nsamp = moldr.util.nsamp_init(nsamp_par, ntaudof)

    run_tau(
        zma=zma,
        spc_info=spc_info,
        thy_level=thy_level,
        nsamp=nsamp,
        tors_range_dct=tors_range_dct,
        tau_run_fs=tau_run_fs,
        tau_save_fs=tau_save_fs,
        script_str=script_str,
        overwrite=overwrite,
        **opt_kwargs,
    )

    save_tau(
        tau_run_fs=tau_run_fs,
        tau_save_fs=tau_save_fs,
    )


def run_tau(
        zma, spc_info, thy_level, nsamp, tors_range_dct,
        tau_run_fs, tau_save_fs, script_str, overwrite, **kwargs):
    """ run sampling algorithm to find tau dependent geometries
    """
    if not tors_range_dct:
        print("No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    tau_save_fs.trunk.create()

    vma = automol.zmatrix.var_(zma)
    if tau_save_fs.trunk.file.vmatrix.exists():
        existing_vma = tau_save_fs.trunk.file.vmatrix.read()
        assert vma == existing_vma
    tau_save_fs.trunk.file.vmatrix.write(vma)
    idx = 0
    nsamp0 = nsamp
    inf_obj = autofile.system.info.tau_trunk(0, tors_range_dct)
    while True:
        if tau_save_fs.trunk.file.info.exists():
            inf_obj_s = tau_save_fs.trunk.file.info.read()
            nsampd = inf_obj_s.nsamp
        elif tau_save_fs.trunk.file.info.exists():
            inf_obj_r = tau_save_fs.trunk.file.info.read()
            nsampd = inf_obj_r.nsamp
        else:
            nsampd = 0

        nsamp = nsamp0 - nsampd
        if nsamp <= 0:
            print('Reached requested number of samples. '
                  'Tau sampling complete.')
            break
        else:
            print("    New nsamp is {:d}.".format(nsamp))

            samp_zma, = automol.zmatrix.samples(zma, 1, tors_range_dct)
            tid = autofile.system.generate_new_tau_id()
            locs = [tid]

            tau_run_fs.leaf.create(locs)
            run_path = tau_run_fs.leaf.path(locs)

            idx += 1
            print("Run {}/{}".format(idx, nsamp0))
            run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=tau_run_fs,
                geom=samp_zma,
                spc_info=spc_info,
                thy_level=thy_level,
                overwrite=overwrite,
                frozen_coordinates=tors_range_dct.keys(),
                **kwargs
            )

            nsampd += 1
            inf_obj.nsamp = nsampd
            tau_save_fs.trunk.file.info.write(inf_obj)
            tau_run_fs.trunk.file.info.write(inf_obj)


def save_tau(tau_run_fs, tau_save_fs):
    """ save the tau dependent geometries that have been found so far
    """

    saved_geos = [tau_save_fs.leaf.file.geometry.read(locs)
                  for locs in tau_save_fs.leaf.existing()]

    if not tau_run_fs.trunk.exists():
        print("No tau geometries to save. Skipping...")
    else:
        for locs in tau_run_fs.leaf.existing():
            run_path = tau_run_fs.leaf.path(locs)
            run_fs = autofile.fs.run(run_path)

            print("Reading from tau run at {}".format(run_path))

            ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)

                geo = elstruct.reader.opt_geometry(prog, out_str)

                save_path = tau_save_fs.leaf.path(locs)
                print(" - Saving...")
                print(" - Save path: {}".format(save_path))

                tau_save_fs.leaf.create(locs)
                tau_save_fs.leaf.file.geometry_info.write(inf_obj, locs)
                tau_save_fs.leaf.file.geometry_input.write(inp_str, locs)
                tau_save_fs.leaf.file.energy.write(ene, locs)
                tau_save_fs.leaf.file.geometry.write(geo, locs)

                saved_geos.append(geo)

        # update the tau trajectory file
        moldr.util.traj_sort(tau_save_fs)


def tau_pf_write(
        name, save_prefix,
        run_grad=False, run_hess=False):
    """ Print out data fle for partition function evaluation 
    """
    cnf_save_fs = autofile.fs.conformer(save_prefix)
    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        ene_ref = cnf_save_fs.leaf.file.energy.read(min_cnf_locs)
        print('ene_ref')
        print(ene_ref)

    tau_save_fs = autofile.fs.tau(save_prefix)
    evr = name+'\n'
    # cycle through saved tau geometries
    idx = 0
    for locs in tau_save_fs.leaf.existing():
        geo = tau_save_fs.leaf.file.geometry.read(locs)
        ene = tau_save_fs.leaf.file.energy.read(locs)
        ene = (ene - ene_ref)*qcc.conversion_factor('hartree', 'kcal/mol')
        ene_str = autofile.file.write.energy(ene)
        geo_str = autofile.file.write.geometry(geo)

        idx += 1
        idx_str = str(idx)

        evr += 'Sampling point'+idx_str+'\n'
        evr += 'Energy'+'\n'
        evr += ene_str+'\n'
        evr += 'Geometry'+'\n'
        evr += geo_str+'\n'
        if run_grad:
            grad = tau_save_fs.leaf.file.gradient.read(locs)
            grad_str = autofile.file.write.gradient(grad)
            evr += 'Gradient'+'\n'
            evr += grad_str
        if run_hess:
            hess = tau_save_fs.leaf.file.hessian.read(locs)
            hess_str = autofile.file.write.hessian(hess)
            evr += 'Hessian'+'\n'
            evr += hess_str+'\n'

    file_name = os.path.join(save_prefix, 'TAU', 'tau.out')
    with open(file_name, 'w') as tau_file:
        tau_file.write(evr)

    temp_list = [300., 500., 750., 1000., 1500.]
    for temp in temp_list:
        sumq = 0.
        sum2 = 0.
        idx = 0
        print('integral convergence for T = ', temp)
        for locs in tau_save_fs.leaf.existing():
            idx += 1
            ene = tau_save_fs.leaf.file.energy.read(locs)
            ene = (ene - ene_ref)*qcc.conversion_factor('hartree', 'kcal/mol')
            tmp = numpy.exp(-ene*349.7/(0.695*temp))
            sumq = sumq + tmp
            sum2 = sum2 + tmp**2
            sigma = numpy.sqrt(
                (abs(sum2/float(idx)-(sumq/float(idx))**2))/float(idx))
            print(sumq/float(idx), sigma, 100.*sigma*float(idx)/sumq, idx)


