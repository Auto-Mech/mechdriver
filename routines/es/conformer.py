""" drivers for conformer
"""

import numpy
import automol
import elstruct
import autofile

# New libs
from routines.es import util
from lib.phydat import phycon
from lib.runner import driver
from lib.runner import par as runpar
from lib.filesystem import minc as fsmin


def conformer_sampling(
        spc_info, thy_level, thy_save_fs, cnf_run_fs, cnf_save_fs, script_str,
        overwrite, saddle=False, nsamp_par=(False, 3, 3, 1, 50, 50),
        tors_names='', dist_info=(), two_stage=False, rxn_class='', **kwargs):
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
        ntaudof = len(
            automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
        nsamp = util.nsamp_init(nsamp_par, ntaudof)
    else:
        ntaudof = len(tors_names)
        nsamp = util.nsamp_init(nsamp_par, ntaudof)

    save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        saddle=saddle,
        dist_info=dist_info,
        rxn_class=rxn_class
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
        dist_info=dist_info,
        rxn_class=rxn_class
    )

    # save information about the minimum energy conformer in top directory
    min_cnf_locs = fsmin.min_energy_conformer_locators(cnf_save_fs)
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


def single_conformer(spc_info, thy_level, filesys, overwrite,
                     saddle=False, dist_info=()):
    """ generate single optimized geometry for
        randomly sampled initial torsional angles
    """
    mc_nsamp = [False, 0, 0, 0, 0, 1]
    sp_script_str, _, kwargs, _ = runpar.run_qchem_par(*thy_level[0:2])
    thy_save_fs = filesys[3]
    two_stage = False
    if saddle:
        two_stage = True
    conformer_sampling(
        spc_info=spc_info,
        thy_level=thy_level,
        thy_save_fs=thy_save_fs,
        cnf_run_fs=filesys[4],
        cnf_save_fs=filesys[5],
        script_str=sp_script_str,
        overwrite=overwrite,
        nsamp_par=mc_nsamp,
        saddle=saddle,
        dist_info=dist_info,
        two_stage=two_stage,
        **kwargs,
    )


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
        print(existing_vma)
        print(vma)
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

            if nsampd > 0:
                samp_zma, = automol.zmatrix.samples(zma, 1, tors_range_dct)
            else:
                samp_zma = zma

            cid = autofile.system.generate_new_conformer_id()
            locs = [cid]

            cnf_run_fs.leaf.create(locs)
            cnf_run_path = cnf_run_fs.leaf.path(locs)
            run_fs = autofile.fs.run(cnf_run_path)

            idx += 1
            print("Run {}/{}".format(nsampd+1, nsamp0))
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
                    thy_level=thy_level,
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
                        thy_level=thy_level,
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
                    thy_level=thy_level,
                    overwrite=overwrite,
                    saddle=saddle,
                    **kwargs
                )

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


def save_conformers(cnf_run_fs, cnf_save_fs, saddle=False,
                    dist_info=(), rxn_class=''):
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
                        print('zma conformer:\n', automol.zmatrix.string(zma))
                        dist_name = dist_info[0]
                        dist_len = dist_info[1]
                        ts_bnd = automol.zmatrix.bond_idxs(zma, dist_name)
                        ts_bnd1 = min(ts_bnd)
                        ts_bnd2 = max(ts_bnd)
                        conf_dist_len = automol.zmatrix.values(zma)[dist_name]
                        brk_name = dist_info[3]
                        cent_atm = None
                        ldist = len(dist_info)
                        print('ldist test:', ldist, dist_info)
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
                            print('angle test in conformer selection:',
                                  angle, conf_ang)
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
                            cls = is_atom_closest_to_bond_atom(
                                zma, ts_bnd2, conf_dist_len)
                            if not cls:
                                print(" - Transition State conformer has",
                                      "diverged from original structure of",
                                      "dist {:.3f} with dist {:.3f}".format(
                                          dist_len, conf_dist_len))
                                print("Radical atom now new nearest neighbor")
                                continue
                            if abs(conf_dist_len - dist_len) > max_disp:
                                print(" - Transition State conformer has",
                                      "diverged from original structure of",
                                      "dist {:.3f} with dist {:.3f}".format(
                                          dist_len, conf_dist_len))
                                continue
                            symbols = automol.zmatrix.symbols(zma)
                            equi_bnd = 0.
                            if symbols[ts_bnd2] == 'H':
                                if symbols[ts_bnd1] == 'H':
                                    equi_bnd = 0.75 * phycon.ANG2BOHR
                                if symbols[ts_bnd1] == 'C':
                                    equi_bnd = 1.09 * phycon.ANG2BOHR
                                elif symbols[ts_bnd1] == 'N':
                                    equi_bnd = 1.01
                                elif symbols[ts_bnd1] == 'O':
                                    equi_bnd = 0.96 * phycon.ANG2BOHR
                            if symbols[ts_bnd2] == 'C':
                                if symbols[ts_bnd1] == 'H':
                                    equi_bnd = 1.09 * phycon.ANG2BOHR
                                if symbols[ts_bnd1] == 'C':
                                    equi_bnd = 1.5 * phycon.ANG2BOHR
                                elif symbols[ts_bnd1] == 'N':
                                    equi_bnd = 1.45 * phycon.ANG2BOHR
                                elif symbols[ts_bnd1] == 'O':
                                    equi_bnd = 1.4 * phycon.ANG2BOHR
                            if symbols[ts_bnd2] == 'N':
                                if symbols[ts_bnd1] == 'H':
                                    equi_bnd = 1.01 * phycon.ANG2BOHR
                                if symbols[ts_bnd1] == 'C':
                                    equi_bnd = 1.45 * phycon.ANG2BOHR
                                elif symbols[ts_bnd1] == 'N':
                                    equi_bnd = 1.4 * phycon.ANG2BOHR
                                elif symbols[ts_bnd1] == 'O':
                                    equi_bnd = 1.35 * phycon.ANG2BOHR
                            if symbols[ts_bnd2] == 'O':
                                if symbols[ts_bnd1] == 'H':
                                    equi_bnd = 0.96 * phycon.ANG2BOHR
                                if symbols[ts_bnd1] == 'C':
                                    equi_bnd = 1.4 * phycon.ANG2BOHR
                                elif symbols[ts_bnd1] == 'N':
                                    equi_bnd = 1.35 * phycon.ANG2BOHR
                                elif symbols[ts_bnd1] == 'O':
                                    equi_bnd = 1.3 * phycon.ANG2BOHR
                            displace_from_equi = conf_dist_len - equi_bnd
                            print('distance_from_equi test:',
                                  conf_dist_len, equi_bnd, dist_len)
                            print('bnd atoms:', ts_bnd1, ts_bnd2,
                                  symbols[ts_bnd1], symbols[ts_bnd2])
                            print('symbols:', symbols[ts_bnd1],
                                  symbols[ts_bnd2])
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
                        if cnf_save_fs.trunk.file.vmatrix.exists():
                            exist_vma = cnf_save_fs.trunk.file.vmatrix.read()
                            if vma != exist_vma:
                                print(" - Isomer is not the same as starting",
                                      "isomer. Skipping...")
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
            if distance < bond_dist:
                atom_closest = False
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
    print('geom: \n', automol.geom.string(geo))
    if automol.geom.is_atom(geo):
        int_sym_num = 1.
    else:
        if not saddle:
            tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
        if tors_names is None:
            int_sym_num = 1.
        else:
            ethrsh = 1.e-5
            locs_lst = cnf_save_fs.leaf.existing()
            int_sym_num = 1.
            if locs_lst:
                enes = [cnf_save_fs.leaf.file.energy.read(locs)
                        for locs in locs_lst]
                geos = [cnf_save_fs.leaf.file.geometry.read(locs)
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
                    print('geo_conf: \n', automol.geom.string(geo_sim_i))
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
                                elif are_torsions_same(new_geo, geo_sim_j):
                                    new_geom = False
                                    break
                        if new_geom:
                            geo_sim2.append(new_geo)
                            int_sym_num += 1
                            print('int_sym_num new:', int_sym_num)
                            print('new_geo: \n', automol.geom.string(new_geo))
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
    print('external/internal test:', ext_sym, int_sym)
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
