""" drivers for conformer
"""

import numpy
import automol
import elstruct
import autofile

# New libs
from routines.es import util
from lib.runner import driver
from lib.runner import par as runpar
from lib.filesystem import minc as fsmin
from lib.phydat import bnd


def conformer_sampling(zma, spc_info,
                       mod_thy_info, thy_save_fs,
                       cnf_run_fs, cnf_save_fs,
                       script_str, overwrite,
                       saddle=False, nsamp_par=(False, 3, 3, 1, 50, 50),
                       tors_names='', dist_info=(),
                       two_stage=False, rxn_class='', **kwargs):
    """ Find the minimum energy conformer by optimizing from nsamp random
    initial torsional states
    """

    ich = spc_info[0]
    coo_names = []

    # Read the geometry and zma from the ini file system
    # if not saddle:
    #     geo = thy_save_fs[-1].file.geometry.read(mod_ini_thy_info[1:4])
    #     tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    #     zma = automol.geom.zmatrix(geo)
    # else:
    #     geo = thy_save_fs[0].file.geometry.read()
    #     zma = thy_save_fs[0].file.zmatrix.read()
    #     coo_names.append(tors_names)
    if saddle:
        coo_names.append(tors_names)

    tors_ranges = tuple((0, 2*numpy.pi) for tors in tors_names)
    tors_range_dct = dict(zip(tors_names, tors_ranges))
    if not saddle:
        gra = automol.inchi.graph(ich)
        ntaudof = len(
            automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
        nsamp = util.nsamp_init(nsamp_par, ntaudof)
    else:
        ntaudof = len(tors_names)
        nsamp = util.nsamp_init(nsamp_par, ntaudof)

    print('\nSaving any conformers in run filesys...')
    save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        thy_info=mod_thy_info,
        saddle=saddle,
        dist_info=dist_info,
        rxn_class=rxn_class
    )

    print('\nSampling for more conformers if needed...')
    run_conformers(
        zma=zma,
        spc_info=spc_info,
        thy_info=mod_thy_info,
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

    print('\nSaving any newly found conformers in run filesys...')
    save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        thy_info=mod_thy_info,
        saddle=saddle,
        dist_info=dist_info,
        rxn_class=rxn_class
    )

    # Save information about the minimum energy conformer in top directory
    min_cnf_locs = fsmin.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        geo = cnf_save_fs[-1].file.geometry.read(min_cnf_locs)
        zma = cnf_save_fs[-1].file.zmatrix.read(min_cnf_locs)
        if not saddle:
            assert automol.zmatrix.almost_equal(zma, automol.geom.zmatrix(geo))
            thy_save_fs[-1].file.geometry.write(geo, mod_thy_info[1:4])
            thy_save_fs[-1].file.zmatrix.write(zma, mod_thy_info[1:4])

        else:
            thy_save_fs[0].file.geometry.write(geo)
            thy_save_fs[0].file.zmatrix.write(zma)


def single_conformer(zma, spc_info, thy_info,
                     thy_save_fs, cnf_run_fs, cnf_save_fs,
                     overwrite, saddle=False, dist_info=()):
    """ generate single optimized geometry for
        randomly sampled initial torsional angles
    """
    opt_script_str, _, kwargs, _ = runpar.run_qchem_par(*thy_info[0:2])
    conformer_sampling(
        zma=zma,
        spc_info=spc_info,
        mod_thy_info=thy_info,
        thy_save_fs=thy_save_fs,
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        script_str=opt_script_str,
        overwrite=overwrite,
        nsamp_par=[False, 0, 0, 0, 0, 1],
        saddle=saddle,
        dist_info=dist_info,
        two_stage=saddle,
        **kwargs,
    )


def run_conformers(
        zma, spc_info, thy_info, nsamp, tors_range_dct,
        cnf_run_fs, cnf_save_fs, script_str, overwrite, saddle, two_stage,
        **kwargs):
    """ run sampling algorithm to find conformers
    """
    if not tors_range_dct:
        print("No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    cnf_save_fs[0].create()
    vma = automol.zmatrix.var_(zma)
    if cnf_save_fs[0].file.vmatrix.exists():
        existing_vma = cnf_save_fs[0].file.vmatrix.read()
        assert vma == existing_vma
    cnf_save_fs[0].file.vmatrix.write(vma)
    idx = 0
    nsamp0 = nsamp
    inf_obj = autofile.system.info.conformer_trunk(0, tors_range_dct)
    if cnf_save_fs[0].file.info.exists():
        inf_obj_s = cnf_save_fs[0].file.info.read()
        nsampd = inf_obj_s.nsamp
    elif cnf_run_fs[0].file.info.exists():
        inf_obj_r = cnf_run_fs[0].file.info.read()
        nsampd = inf_obj_r.nsamp
    else:
        nsampd = 0
   
    tot_samp = nsamp - nsampd
    print('Number of samples requested:', nsamp)
    print('Number of samples that have been currently run:', nsampd, '\n')

    samp_idx = 1
    while True:
        nsamp = nsamp0 - nsampd
        # Break the while loop if enough sampls completed
        if nsamp <= 0:
            print('Reached requested number of samples. '
                  'Conformer search complete.')
            break

        # Run the conformer sampling
        if nsampd > 0:
            samp_zma, = automol.zmatrix.samples(zma, 1, tors_range_dct)
        else:
            samp_zma = zma

        cid = autofile.system.generate_new_conformer_id()
        locs = [cid]

        cnf_run_fs[-1].create(locs)
        cnf_run_path = cnf_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(cnf_run_path)

        print("Run {}/{}".format(samp_idx, tot_samp))
        tors_names = list(tors_range_dct.keys())
        if two_stage and tors_names:
            print('Stage one beginning, holding the coordinates constant',
                  tors_names)
            driver.run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geom=samp_zma,
                spc_info=spc_info,
                thy_level=thy_info,
                overwrite=overwrite,
                frozen_coordinates=[tors_names],
                saddle=saddle,
                **kwargs
            )
            print('Stage one success, reading for stage 2')
            ret = driver.read_job(
                job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret:
                sinf_obj, _, out_str = ret
                prog = sinf_obj.prog
                samp_zma = elstruct.reader.opt_zmatrix(prog, out_str)
                print('Stage one success beginning stage two on', samp_zma)
                driver.run_job(
                    job=elstruct.Job.OPTIMIZATION,
                    script_str=script_str,
                    run_fs=run_fs,
                    geom=samp_zma,
                    spc_info=spc_info,
                    thy_level=thy_info,
                    overwrite=overwrite,
                    saddle=saddle,
                    **kwargs
                )
        else:
            driver.run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geom=samp_zma,
                spc_info=spc_info,
                thy_level=thy_info,
                overwrite=overwrite,
                saddle=saddle,
                **kwargs
            )

        if cnf_save_fs[0].file.info.exists():
            inf_obj_s = cnf_save_fs[0].file.info.read()
            nsampd = inf_obj_s.nsamp
        elif cnf_run_fs[0].file.info.exists():
            inf_obj_r = cnf_run_fs[0].file.info.read()
            nsampd = inf_obj_r.nsamp
        nsampd += 1
        samp_idx += 1
        inf_obj.nsamp = nsampd
        cnf_save_fs[0].file.info.write(inf_obj)
        cnf_run_fs[0].file.info.write(inf_obj)


def save_conformers(cnf_run_fs, cnf_save_fs, thy_info, saddle=False,
                    dist_info=(), rxn_class=''):
    """ save the conformers that have been found so far
    """

    locs_lst = cnf_save_fs[-1].existing()
    seen_geos = [cnf_save_fs[-1].file.geometry.read(locs)
                 for locs in locs_lst]
    seen_enes = [cnf_save_fs[-1].file.energy.read(locs)
                 for locs in locs_lst]

    if not cnf_run_fs[0].exists():
        print("No conformers in run filesys to save.")
    else:
        print("Found conformers in run filesys at level of theory to save.")
        for locs in cnf_run_fs[-1].existing():
            # # Only go through save procedure if conf not in save
            # # may need to get geo, ene, etc; maybe make function
            # if cnf_save_fs[-1].exists(locs):
            #     continue
            # else:
            #     print('New conformer to save...')
            cnf_run_path = cnf_run_fs[-1].path(locs)
            run_fs = autofile.fs.run(cnf_run_path)
            print("Reading from conformer run at {}".format(cnf_run_path))

            ret = driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                if not saddle:
                    gra = automol.geom.graph(geo)
                    conns = automol.graph.connected_components(gra)
                    lconns = len(conns)
                else:
                    lconns = 1
                if lconns > 1:
                    print(" - Geometry is disconnected.. Skipping...")
                else:
                    if saddle:
                        zma = elstruct.reader.opt_zmatrix(prog, out_str)
                        dist_name = dist_info[0]
                        dist_len = dist_info[1]
                        ts_bnd = automol.zmatrix.bond_idxs(zma, dist_name)
                        ts_bnd1 = min(ts_bnd)
                        ts_bnd2 = max(ts_bnd)
                        conf_dist_len = automol.zmatrix.values(zma)[dist_name]
                        brk_name = dist_info[3]
                        cent_atm = None
                        ldist = len(dist_info)
                        if dist_name and brk_name and ldist > 4:
                            angle = dist_info[4]
                            brk_bnd = automol.zmatrix.bond_idxs(zma, brk_name)
                            ang_atms = [0, 0, 0]
                            cent_atm = list(set(brk_bnd) & set(ts_bnd))
                            if cent_atm:
                                ang_atms[1] = cent_atm[0]
                                for idx in brk_bnd:
                                    if idx != ang_atms[1]:
                                        ang_atms[0] = idx
                                for idx in ts_bnd:
                                    if idx != ang_atms[1]:
                                        ang_atms[2] = idx
                                geom = automol.zmatrix.geometry(zma)
                                conf_ang = automol.geom.central_angle(
                                    geom, *ang_atms)
                        max_disp = 0.6
                        if 'addition' in rxn_class:
                            max_disp = 0.8
                        if 'abstraction' in rxn_class:
                            max_disp = 1.4

                        # check forming bond angle similar to ini config
                        if cent_atm and 'elimination' not in rxn_class:
                            # print('angle test in conformer selection:',
                            #       angle, conf_ang)
                            if abs(conf_ang - angle) > .44:
                                print(" - Transition State conformer has",
                                      "diverged from original structure of",
                                      "angle {:.3f} with angle {:.3f}".format(
                                          angle, conf_ang))
                                continue
                        # check if radical atom is closer to some atom
                        # other than the bonding atom
                        if 'add' in rxn_class or 'abst' in rxn_class:
                            print('it is an addition or an abstraction:')
                            print('ts_bnd2 test:', ts_bnd2, zma)
                            cls = is_atom_closest_to_bond_atom(
                                zma, ts_bnd2, conf_dist_len)
                            if not cls:
                                print(" - Transition State conformer has",
                                      "diverged from original structure of",
                                      "dist {:.3f} with dist {:.3f}".format(
                                          dist_len, conf_dist_len))
                                print("Radical atom now has a new nearest neighbor")
                                continue
                            if abs(conf_dist_len - dist_len) > max_disp:
                                print(" - Transition State conformer has",
                                      "diverged from original structure of",
                                      "dist {:.3f} with dist {:.3f}".format(
                                          dist_len, conf_dist_len))
                                continue
                            symbols = automol.zmatrix.symbols(zma)

                            # Set standard equivalent bond len for rxn coord
                            symbols = automol.zmatrix.symbols(zma)
                            symb1, symb2 = symbols[ts_bnd1], symbols[ts_bnd2]
                            if (symb1, symb2) in bnd.LEN_DCT:
                                equi_bnd = bnd.LEN_DCT[(symb1, symb2)]
                            elif (symb2, symb1) in bnd.LEN_DCT:
                                equi_bnd = bnd.LEN_DCT[(symb2, symb1)]
                            else:
                                equi_bnd = 0.0
                            displace_from_equi = conf_dist_len - equi_bnd
                            dchk1 = abs(conf_dist_len - dist_len) > 0.2
                            dchk2 = displace_from_equi < 0.2
                            if dchk1 and dchk2:
                                print(" - Transition State conformer has",
                                      "converged to an",
                                      "equilibrium structure with dist",
                                      " {:.3f} comp with equil {:.3f}".format(
                                          conf_dist_len, equi_bnd))
                                continue
                        else:
                            if abs(conf_dist_len - dist_len) > 0.4:
                                print(" - Transition State conformer has",
                                      "diverged from original structure of",
                                      "dist {:.3f} with dist {:.3f}".format(
                                          dist_len, conf_dist_len))
                                continue
                    else:
                        zma = automol.geom.zmatrix(geo)
                    unique = is_unique_tors_dist_mat_energy(
                        geo, ene, seen_geos, seen_enes, saddle)

                    if not unique:
                        print(" - Geometry is not unique. Skipping...")
                    else:
                        vma = automol.zmatrix.var_(zma)
                        if cnf_save_fs[0].file.vmatrix.exists():
                            exist_vma = cnf_save_fs[0].file.vmatrix.read()
                            if vma != exist_vma:
                                print(" - Isomer is not the same as starting",
                                      "isomer. Skipping...")
                            else:
                                save_path = cnf_save_fs[-1].path(locs)
                                print(" - Geometry is unique. Saving...")
                                print(" - Save path: {}".format(save_path))

                                cnf_save_fs[-1].create(locs)
                                cnf_save_fs[-1].file.geometry_info.write(
                                    inf_obj, locs)
                                cnf_save_fs[-1].file.geometry_input.write(
                                    inp_str, locs)
                                cnf_save_fs[-1].file.energy.write(ene, locs)
                                cnf_save_fs[-1].file.geometry.write(geo, locs)
                                cnf_save_fs[-1].file.zmatrix.write(zma, locs)

                                # Saving the energy to am SP filesys
                                print(" - Saving energy...")
                                sp_save_fs = autofile.fs.single_point(
                                    save_path)
                                sp_save_fs[-1].create(thy_info[1:4])
                                sp_save_fs[-1].file.input.write(
                                    inp_str, thy_info[1:4])
                                sp_save_fs[-1].file.info.write(
                                    inf_obj, thy_info[1:4])
                                sp_save_fs[-1].file.energy.write(
                                    ene, thy_info[1:4])

                    seen_geos.append(geo)
                    seen_enes.append(ene)

        # Update the conformer trajectory file
        print('')
        fsmin.traj_sort(cnf_save_fs)


def is_atom_closest_to_bond_atom(zma, idx_rad, bond_dist):
    """ Check to see whether the radical atom is still closest to the bond
        formation site.
    """
    geo = automol.zmatrix.geometry(zma)
    atom_closest = True
    for idx, _ in enumerate(geo):
        if idx < idx_rad:
            distance = automol.geom.distance(geo, idx, idx_rad)
            if distance < bond_dist-0.01:
                atom_closest = False
                print('idx test:', idx, distance, bond_dist)
    return atom_closest


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


def int_sym_num_from_sampling(
        geo, ene, cnf_save_fs, saddle=False, frm_bnd_key=(),
        brk_bnd_key=(), form_coords=(), tors_names=()):
    """ Determine the symmetry number for a given conformer geometry.
    (1) Explore the saved conformers to find the list of similar conformers -
        i.e. those with a coulomb matrix and energy that are equivalent
        to those for the reference geometry.
    (2) Expand each of those similar conformers by applying
        rotational permutations to each of the terminal groups.
    (3) Count how many distinct distance matrices there are in
        the fully expanded conformer list.
    """

    # Note: ignoring for saddle points the possibility that two configurations
    # differ only in their torsional values.
    # As a result, the symmetry factor is a lower bound of the true value
    if automol.geom.is_atom(geo):
        int_sym_num = 1.
    else:
        if not saddle:
            tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
        if tors_names is None:
            int_sym_num = 1.
        else:
            ethrsh = 1.e-5
            locs_lst = cnf_save_fs[-1].existing()
            int_sym_num = 1.
            if locs_lst:
                enes = [cnf_save_fs[-1].file.energy.read(locs)
                        for locs in locs_lst]
                geos = [cnf_save_fs[-1].file.geometry.read(locs)
                        for locs in locs_lst]
                geo_sim = []
                geo_sim2 = []
                ene_sim = []
                for geoi, enei in zip(geos, enes):
                    if enei - enes[0] < ethrsh:
                        geo_lst = [geoi]
                        ene_lst = [enei]
                        unique = is_unique_coulomb_energy(
                            geo, ene, geo_lst, ene_lst)
                        if not unique:
                            geo_sim.append(geoi)
                            ene_sim.append(enei)

                int_sym_num = 0
                for geo_sim_i in geo_sim:
                    new_geos = automol.geom.rot_permutated_geoms(
                        geo_sim_i, saddle,
                        frm_bnd_key, brk_bnd_key, form_coords)
                    for new_geo in new_geos:
                        new_geom = True
                        for geo_sim_j in geo_sim2:
                            if automol.geom.almost_equal_dist_mat(
                                    new_geo, geo_sim_j, thresh=3e-1):
                                if saddle:
                                    new_geom = False
                                    break
                                if are_torsions_same(new_geo, geo_sim_j):
                                    new_geom = False
                                    break
                        if new_geom:
                            geo_sim2.append(new_geo)
                            int_sym_num += 1
    return int_sym_num


def symmetry_factor(
        geo, ene, cnf_save_fs, saddle=False, frm_bnd_key=(), brk_bnd_key=(),
        form_coords=(), tors_names=()):
    """ obtain overall symmetry factor for a geometry as a product
        of the external symmetry factor and the internal symmetry number
    """
    # Note: ignoring for saddle points the possibility that two configurations
    # differ only in their torsional values.
    # As a result, the symmetry factor is a lower bound of the true value
    ext_sym = automol.geom.external_symmetry_factor(geo)
    if not saddle:
        tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    if tors_names:
        int_sym = int_sym_num_from_sampling(
            geo, ene, cnf_save_fs, saddle,
            frm_bnd_key, brk_bnd_key,
            form_coords, tors_names)
    else:
        int_sym = 1
    sym_fac = ext_sym * int_sym
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
    dtol = 0.09
    same_dihed = True
    zma = automol.geom.zmatrix(geo)
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    zmai = automol.geom.zmatrix(geoi)
    tors_namesi = automol.geom.zmatrix_torsion_coordinate_names(geoi)
    for idx, tors_name in enumerate(tors_names):
        val = automol.zmatrix.values(zma)[tors_name]
        vali = automol.zmatrix.values(zmai)[tors_namesi[idx]]
        valip = vali+2.*numpy.pi
        valim = vali-2.*numpy.pi
        vchk1 = abs(val - vali)
        vchk2 = abs(val - valip)
        vchk3 = abs(val - valim)
        if vchk1 > dtol and vchk2 > dtol and vchk3 > dtol:
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
                    geo, geoi, thresh=3e-1):
                # check dihedrals
                if saddle:
                    unique = False
                elif are_torsions_same(geo, geoi):
                    unique = False
    return unique
