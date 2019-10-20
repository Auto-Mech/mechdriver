""" drivers for conformer
"""

import numpy
from qcelemental import constants as qcc
import automol
import elstruct
import autofile
import moldr
from automol.convert import _pyx2z

WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

def conformer_sampling(
        spc_info, thy_level, thy_save_fs, cnf_run_fs, cnf_save_fs, script_str,
        overwrite, saddle=False, nsamp_par=(False, 3, 3, 1, 50, 50),
        tors_names='', dist_info=[], two_stage=False, **kwargs):
    """ Find the minimum energy conformer by optimizing from nsamp random
    initial torsional states
    """
    ich = spc_info[0]
    coo_names = []
    if not saddle:
        geo = thy_save_fs.leaf.file.geometry.read(thy_level[1:4])
        tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
        zma = automol.geom.zmatrix(geo)
    else:
        geo = thy_save_fs.trunk.file.geometry.read()
        zma = thy_save_fs.trunk.file.zmatrix.read()
        coo_names.append(tors_names)


    tors_ranges = tuple((0, 2*numpy.pi) for tors in tors_names)
    tors_range_dct = dict(zip(tors_names, tors_ranges))
    if not saddle:
        gra = automol.inchi.graph(ich)
        ntaudof = len(automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
        nsamp = moldr.util.nsamp_init(nsamp_par, ntaudof)
    else:
        ntaudof = len(tors_names)
        nsamp = moldr.util.nsamp_init(nsamp_par, ntaudof)
    save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        saddle=saddle,
        dist_info=dist_info
    )

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
        saddle=saddle,
        two_stage=two_stage,
        **kwargs,
    )
    save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        saddle=saddle,
        dist_info=dist_info
    )

    # save information about the minimum energy conformer in top directory
    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        zma = cnf_save_fs.leaf.file.zmatrix.read(min_cnf_locs)
        if not saddle:
            assert automol.zmatrix.almost_equal(zma, automol.geom.zmatrix(geo))
            thy_save_fs.leaf.file.geometry.write(geo, thy_level[1:4])
            thy_save_fs.leaf.file.zmatrix.write(zma, thy_level[1:4])

        else:
            thy_save_fs.trunk.file.geometry.write(geo)
            thy_save_fs.trunk.file.zmatrix.write(zma)

        ene = cnf_save_fs.leaf.file.energy.read(min_cnf_locs)

def run_conformers(
        zma, spc_info, thy_level, nsamp, tors_range_dct,
        cnf_run_fs, cnf_save_fs, script_str, overwrite, saddle, two_stage,
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
    if cnf_save_fs.trunk.file.info.exists():
        inf_obj_s = cnf_save_fs.trunk.file.info.read()
        nsampd = inf_obj_s.nsamp
    elif cnf_run_fs.trunk.file.info.exists():
        inf_obj_r = cnf_run_fs.trunk.file.info.read()
        nsampd = inf_obj_r.nsamp
    else:
        nsampd = 0

    while True:
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

            idx += 1
            print("Run {}/{}".format(nsampd+1, nsamp0))
            tors_names = list(tors_range_dct.keys())
            if two_stage and len(tors_names) > 0:
                print('Stage one beginning, holding the coordinates constant', tors_names, samp_zma)
                moldr.driver.run_job(
                    job=elstruct.Job.OPTIMIZATION,
                    script_str=script_str,
                    run_fs=run_fs,
                    geom=samp_zma,
                    spc_info=spc_info,
                    thy_level=thy_level,
                    overwrite=overwrite,
                    frozen_coordinates=[tors_names],
                    saddle=saddle,
                    **kwargs
                )
                print('Stage one success, reading for stage 2')
                ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
                if ret:
                    sinf_obj, inp_str, out_str = ret
                    prog = sinf_obj.prog
                    samp_zma = elstruct.reader.opt_zmatrix(prog, out_str)
                    print('Stage one success beginning stage two on', samp_zma)
                    moldr.driver.run_job(
                        job=elstruct.Job.OPTIMIZATION,
                        script_str=script_str,
                        run_fs=run_fs,
                        geom=samp_zma,
                        spc_info=spc_info,
                        thy_level=thy_level,
                        overwrite=overwrite,
                        saddle=saddle,
                        **kwargs
                    )
            else:
                moldr.driver.run_job(
                    job=elstruct.Job.OPTIMIZATION,
                    script_str=script_str,
                    run_fs=run_fs,
                    geom=samp_zma,
                    spc_info=spc_info,
                    thy_level=thy_level,
                    overwrite=overwrite,
                    saddle=saddle,
                    **kwargs
                )
            #moldr.driver.run_job(
            #    job=elstruct.Job.OPTIMIZATION,
            #    script_str=script_str,
            #    run_fs=run_fs,
            #    geom=samp_zma,
            #    spc_info=spc_info,
            #    thy_level=thy_level,
            #    overwrite=overwrite,
            #    **kwargs
            #)

            if cnf_save_fs.trunk.file.info.exists():
                inf_obj_s = cnf_save_fs.trunk.file.info.read()
                nsampd = inf_obj_s.nsamp
            elif cnf_run_fs.trunk.file.info.exists():
                inf_obj_r = cnf_run_fs.trunk.file.info.read()
                nsampd = inf_obj_r.nsamp
            nsampd += 1
            inf_obj.nsamp = nsampd
            cnf_save_fs.trunk.file.info.write(inf_obj)
            cnf_run_fs.trunk.file.info.write(inf_obj)


def save_conformers(cnf_run_fs, cnf_save_fs, saddle=False, dist_info=[]):
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
                if saddle:
                    gra = automol.geom.weakly_connected_graph(geo)
                else:
                    gra = automol.geom.graph(geo)

                conns = automol.graph.connected_components(gra)
                if len(conns) >  1:
                    print(" - Geometry is disconnected.. Skipping...")
                else:
                    if saddle:
                        zma = elstruct.reader.opt_zmatrix(prog, out_str)
                        dist_name = dist_info[0]
                        dist_len = dist_info[1]
                        conf_dist_len = automol.zmatrix.values(zma)[dist_name]
                        print('distance test:', dist_len, conf_dist_len)
                        if abs(conf_dist_len - dist_len) > 0.3:
                            print(" - Transition State conformer has diverged from original structure of dist {:.3f} with dist {:.3f}".format(dist_len, conf_dist_len))
                            continue
                    else:
                        zma = automol.geom.zmatrix(geo)
                    unique = is_unique_tors_dist_mat_energy(geo, ene, seen_geos, seen_enes, saddle)

                    if not unique:
                        print(" - Geometry is not unique. Skipping...")
                    else:
                        vma = automol.zmatrix.var_(zma)
                        if cnf_save_fs.trunk.file.vmatrix.exists():
                            existing_vma = cnf_save_fs.trunk.file.vmatrix.read()
                            if vma != existing_vma:
                                print(" - Isomer is not the same as starting isomer. Skipping...")
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
        moldr.util.traj_sort(cnf_save_fs)


def run_gradient(
        spc_info, thy_level, geo_run_fs, geo_save_fs, locs,
        script_str, overwrite, **kwargs):
    """ Determine the gradient for the geometry in the given location
    """

    geo_run_path = geo_run_fs.leaf.path(locs)
    geo_save_path = geo_save_fs.leaf.path(locs)
    geo = geo_save_fs.leaf.file.geometry.read(locs)
    run_fs = autofile.fs.run(geo_run_path)

    print('Running gradient')
    moldr.driver.run_job(
        job='gradient',
        script_str=script_str,
        run_fs=run_fs,
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

        if automol.geom.is_atom(geo):
            freqs = ()
        else:
            print(" - Reading gradient from output...")
            grad = elstruct.reader.gradient(inf_obj.prog, out_str)

            print(" - Saving gradient...")
            print(" - Save path: {}".format(geo_save_path))
            geo_save_fs.leaf.file.gradient_info.write(inf_obj, locs)
            geo_save_fs.leaf.file.gradient_input.write(inp_str, locs)
            geo_save_fs.leaf.file.gradient.write(grad, locs)


def run_hessian(
        spc_info, thy_level, geo_run_fs, geo_save_fs, locs,
        script_str, overwrite, **kwargs):
    """ Determine the hessian for the geometry in the given location
    """

    geo_run_path = geo_run_fs.leaf.path(locs)
    geo_save_path = geo_save_fs.leaf.path(locs)
    geo = geo_save_fs.leaf.file.geometry.read(locs)
    run_fs = autofile.fs.run(geo_run_path)

    print('Running hessian')
    moldr.driver.run_job(
        job='hessian',
        script_str=script_str,
        run_fs=run_fs,
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

            print(" - Saving hessian...")
            print(" - Save path: {}".format(geo_save_path))
            geo_save_fs.leaf.file.hessian_info.write(inf_obj, locs)
            geo_save_fs.leaf.file.hessian_input.write(inp_str, locs)
            geo_save_fs.leaf.file.hessian.write(hess, locs)
            geo_save_fs.leaf.file.harmonic_frequencies.write(freqs, locs)


def run_vpt2(
        spc_info, thy_level, geo_run_fs, geo_save_fs, locs,
        script_str, overwrite, **kwargs):
    """ Perform vpt2 analysis for the geometry in the given location
    """

    geo_run_path = geo_run_fs.leaf.path(locs)
    geo_save_path = geo_save_fs.leaf.path(locs)
    geo = geo_save_fs.leaf.file.geometry.read(locs)
    run_fs = autofile.fs.run(geo_run_path)

    print('Running vpt2')
    moldr.driver.run_job(
        job='vpt2',
        script_str=script_str,
        run_fs=run_fs,
        geom=geo,
        spc_info=spc_info,
        thy_level=thy_level,
        overwrite=overwrite,
        **kwargs,
    )

    ret = moldr.driver.read_job(
        job='vpt2',
        run_fs=run_fs,
    )

    if ret is not None:
        inf_obj, inp_str, out_str = ret

        if not automol.geom.is_atom(geo):
            print(" - Reading vpt2 from output...")
            vpt2 = elstruct.reader.vpt2(inf_obj.prog, out_str)

            print(" - Saving vpt2...")
            print(" - Save path: {}".format(geo_save_path))
            geo_save_fs.leaf.file.vpt2.info.write(inf_obj, locs)
            geo_save_fs.leaf.file.vpt2.input.write(inp_str, locs)
            geo_save_fs.leaf.file.vpt2.write(vpt2, locs)


def is_unique_coulomb_energy(geo, ene, geo_list, ene_list):
    """ compare given geo with list of geos all to see if any have the same
    coulomb spectrum and energy
    """
    unique = True
    for idx, geoi in enumerate(geo_list):
        enei = ene_list[idx]
        etol = 2.e-5
        if abs(ene-enei) < etol:
            if automol.geom.almost_equal_coulomb_spectrum(
                    geo, geoi, rtol=1e-2):
                unique = False
    return unique


def is_unique_dist_mat_energy(geo, ene, geo_list, ene_list):
    """ compare given geo with list of geos all to see if any have the same
    distance matrix and energy
    """
    unique = True
    for idx, geoi in enumerate(geo_list):
        enei = ene_list[idx]
        etol = 2.e-5
        if abs(ene-enei) < etol:
            if automol.geom.almost_equal_dist_mat(
                    geo, geoi, thresh=1e-1):
                unique = False
    return unique


def int_sym_num_from_sampling(geo, ene, cnf_save_fs, saddle=False, form_coords=[], tors_names=[]):
    """ Determine the symmetry number for a given conformer geometry.
    First explore the saved conformers to find the list of similar conformers -
    i.e. those with a coulomb matrix and energy that are equivalent to those for the
    reference geometry. Then expand each of those similar conformers by applying
    rotational permutations to each of the terminal groups. Finally count how many
    distinct distance matrices there are in the fully expanded conformer list.
    """

    # Note, for now we are ignoring for saddle points the possibility that two configurations
    # differ only in their torsional values. As a result, the symmetry factor is a lower bound on the 
    # true symmetry factor
    if automol.geom.is_atom(geo):
        int_sym_num = 1.
    else:
        if not saddle:
            tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
        else:
            print('tors_names test:', tors_names, len(tors_names))
        if len(tors_names) == 0:
            int_sym_num = 1.
        else:
            ethrsh = 1.e-5
            locs_lst = cnf_save_fs.leaf.existing()
            geo_sim = [geo]
            geo_sim_exp = [geo]
            ene_sim = [ene]
            int_sym_num = 1.
            if locs_lst:
                enes = [cnf_save_fs.leaf.file.energy.read(locs)
                        for locs in locs_lst]
                geos = [cnf_save_fs.leaf.file.geometry.read(locs)
                        for locs in locs_lst]
                for geoi, enei in zip(geos, enes):
                    if enei - enes[0] < ethrsh:
                        geo_lst = [geoi]
                        ene_lst = [enei]
                        if not is_unique_coulomb_energy(geo, ene, geo_lst, ene_lst):
                            geo_sim.append(geoi)
                            ene_sim.append(enei)

                int_sym_num = 0
                for idx_i, geo_sim_i in enumerate(geo_sim):
                    new_geos = automol.geom.rot_permutated_geoms(geo_sim_i, saddle, form_coords)
                    new_geom = True
                    for new_geo in new_geos:
                        #geo_sim_exp.append(new_geo)
                        for idx_j, geo_sim_j in enumerate(geo_sim):
                            if idx_j < idx_i:
                                if automol.geom.almost_equal_dist_mat(
                                        new_geo, geo_sim_j, thresh=1e-1):
                                    if saddle:
                                        new_geom = False 
                                        break
                                    elif are_torsions_same(new_geo, geo_sim_j):
                                        new_geom = False
                                        break
                            else:
                                break
                        if not new_geom:
                            break
                    if new_geom:
                        int_sym_num += len(new_geos)
    return int_sym_num


def external_symmetry_factor(geo):
    """ obtain external symmetry number for a geometry using x2z
    """
    # Get initial external symmetry number
    if automol.geom.is_atom(geo):
        ext_sym_fac = 1.
    else:
        oriented_geom = _pyx2z.to_oriented_geometry(geo)
        ext_sym_fac = oriented_geom.sym_num()
        # print('initial ext_sym_fac:', ext_sym_fac)
        # Change symmetry number if geometry has enantiomers
        if oriented_geom.is_enantiomer():
            ext_sym_fac *= 0.5
        # print('final ext_sym_fac:', ext_sym_fac)
    return ext_sym_fac


def symmetry_factor(geo, ene, cnf_save_fs, saddle=False, form_coords=[], tors_names=[]):
    """ obtain overall symmetry number for a geometry as a product
    of the external and internal symmetry numbers
    """
    # Note, for now we are ignoring for saddle points the possibility that two configurations
    # differ only in their torsional values. As a result, the symmetry factor is a lower bound on the 
    # true symmetry factor
    ext_sym = external_symmetry_factor(geo)
    if not saddle:
        tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    else:
        print('tors_names test:', tors_names, len(tors_names))
    if len(tors_names) > 0:
        int_sym = int_sym_num_from_sampling(geo, ene, cnf_save_fs, saddle, form_coords, tors_names)
    else:
        int_sym = 1
    sym_fac = ext_sym * int_sym
    print('sym test:', ext_sym, int_sym)
    return sym_fac


def is_unique_stereo_dist_mat_energy(geo, ene, geo_list, ene_list):
    """ compare given geo with list of geos all to see if any have the same
    distance matrix and energy and stereo specific inchi
    """
    unique = True
    ich = automol.convert.geom.inchi(geo)
    for idx, geoi in enumerate(geo_list):
        enei = ene_list[idx]
        etol = 2.e-5
        ichi = automol.convert.geom.inchi(geoi)
        # check energy
        if abs(ene-enei) < etol:
            # check distance matrix
            if automol.geom.almost_equal_dist_mat(
                    geo, geoi, thresh=1e-1):
                # check stereo by generates stero label
                ichi = automol.convert.geom.inchi(geoi)
                if ich == ichi:
                    unique = False
    return unique


def are_torsions_same(geo, geoi):
    """ compare all torsional angle values
    """
    dtol = 0.03
    same_dihed = True
    zma = automol.geom.zmatrix(geo)
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    zmai = automol.geom.zmatrix(geoi)
    tors_namesi = automol.geom.zmatrix_torsion_coordinate_names(geoi)
    for idx, tors_name in enumerate(tors_names):
        val = automol.zmatrix.values(zma)[tors_name]
        vali = automol.zmatrix.values(zmai)[tors_namesi[idx]]
        if abs(val - vali) > dtol:
            same_dihed = False
    return same_dihed


def is_unique_tors_dist_mat_energy(geo, ene, geo_list, ene_list, saddle):
    """ compare given geo with list of geos all to see if any have the same
    coulomb spectrum and energy and stereo specific inchi
    """
    unique = True
    etol = 2.e-5
    for idx, geoi in enumerate(geo_list):
        enei = ene_list[idx]
        # check energy
        if abs(ene-enei) < etol:
            # check distance matrix
            if automol.geom.almost_equal_dist_mat(
                    geo, geoi, thresh=1e-1):
                # check dihedrals
                if saddle:
                    unique = False 
                elif are_torsions_same(geo, geoi):
                    unique = False
    return unique
