""" es_runners for conformer
"""

import numpy
import automol
import elstruct
import autofile
from autofile import fs
from phydat import bnd
from mechroutines.es._routines import _util as util
from mechroutines.es import runner as es_runner
from mechlib import filesys
from mechlib.structure import geom as geomprep
from mechlib.structure import ts as tsprep
from mechlib.amech_io import printer as ioprinter


def conformer_sampling(zma, spc_info,
                       mod_thy_info, thy_save_fs,
                       cnf_run_fs, cnf_save_fs,
                       script_str, overwrite,
                       saddle=False, nsamp_par=(False, 3, 3, 1, 50, 50),
                       tors_names='',
                       two_stage=False, retryfail=True,
                       rxn_class='', **kwargs):
    """ Find the minimum energy conformer by optimizing from nsamp random
    initial torsional states
    """

    ich = spc_info[0]
    coo_names = []

    # Read the geometry and zma from the ini file system
    if saddle:
        coo_names.append(tors_names)

    tors_ranges = tuple((0, 2*numpy.pi) for tors in tors_names)
    tors_range_dct = dict(zip(tors_names, tors_ranges))
    ioprinter.debug_message('tors_names', tors_names)
    ioprinter.debug_message('tors_range_dct', tors_range_dct)
    if not saddle:
        gra = automol.inchi.graph(ich)
        ntaudof = len(
            automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
    else:
        ntaudof = len(tors_names)
    nsamp = util.nsamp_init(nsamp_par, ntaudof)

    # Check samples and if nsamp met and no resave

    ioprinter.info_message(
        'Saving any conformers in run filesys...', newline=1)
    save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        thy_info=mod_thy_info,
        saddle=saddle,
        rxn_class=rxn_class,
        orig_ich=ich
    )

    ioprinter.info_message(
        'Sampling for more conformers if needed...', newline=1)
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
        retryfail=retryfail,
        **kwargs,
    )

    ioprinter.info_message(
        'Saving any newly found conformers in run filesys...', newline=1)
    save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        thy_info=mod_thy_info,
        saddle=saddle,
        rxn_class=rxn_class,
        orig_ich=ich
    )

    # Save information about the minimum energy conformer in top directory
    min_cnf_locs, _ = filesys.mincnf.min_energy_conformer_locators(
        cnf_save_fs, mod_thy_info)
    if min_cnf_locs:
        ioprinter.debug_message('min_cnf_locs test in save_conformer:', min_cnf_locs)
        geo = cnf_save_fs[-1].file.geometry.read(min_cnf_locs)
        if not saddle:
            # print(automol.zmatrix.string(zma))
            # print(automol.zmatrix.string(automol.geom.zmatrix(geo)))
            # assert automol.zmatrix.almost_equal(
            #  zma, automol.geom.zmatrix(geo))
            thy_save_fs[-1].file.geometry.write(geo, mod_thy_info[1:4])
            # thy_save_fs[-1].file.zmatrix.write(zma, mod_thy_info[1:4])
        else:
            # thy_save_fs[0].file.zmatrix.write(geo)
            thy_save_fs[0].file.geometry.write(geo)

    return bool(min_cnf_locs)


def single_conformer(zma, spc_info, mod_thy_info,
                     cnf_run_fs, cnf_save_fs,
                     script_str, overwrite,
                     retryfail=True, saddle=False, **kwargs):
    """ generate single optimized geometry to be saved into a
        filesystem
    """

    # Build the filesystem
    locs = [autofile.schema.generate_new_conformer_id()]
    cnf_run_fs[-1].create(locs)
    cnf_run_path = cnf_run_fs[-1].path(locs)
    run_fs = autofile.fs.run(cnf_run_path)

    # Run the optimization
    ioprinter.info_message('Optimizing a single conformer')
    es_runner.run_job(
        job=elstruct.Job.OPTIMIZATION,
        script_str=script_str,
        run_fs=run_fs,
        geom=zma,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        overwrite=overwrite,
        frozen_coordinates=(),
        saddle=saddle,
        retryfail=retryfail,
        **kwargs
    )
    # print('Stage one success, reading for stage 2')
    success, ret = es_runner.read_job(
        job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)

    if success:
        inf_obj, _, out_str = ret
        prog = inf_obj.prog
        method = inf_obj.method
        ene = elstruct.reader.energy(prog, method, out_str)
        geo = elstruct.reader.opt_geometry(prog, out_str)
        zma = elstruct.reader.opt_zmatrix(prog, out_str)
        saved_locs, saved_geos, saved_enes = _saved_cnf_info(
            cnf_save_fs, mod_thy_info)

        if _geo_unique(geo, ene, saved_geos, saved_enes, saddle):
            sym_id = _sym_unique(
                geo, ene, saved_geos, saved_enes)
            if sym_id is None:
                _save_unique_conformer(
                    ret, mod_thy_info, cnf_save_fs, locs,
                    saddle=saddle, zma_locs=(0,))
                saved_geos.append(geo)
                saved_enes.append(ene)
                saved_locs.append(locs)

        # Update the conformer trajectory file
        ioprinter.obj('vspace')
        filesys.mincnf.traj_sort(cnf_save_fs, mod_thy_info)


def run_conformers(
        zma, spc_info, thy_info, nsamp, tors_range_dct,
        cnf_run_fs, cnf_save_fs, script_str, overwrite,
        saddle, two_stage, retryfail,
        **kwargs):
    """ run sampling algorithm to find conformers
    """
    if not tors_range_dct:
        ioprinter.info_message(" - No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    cnf_save_fs[0].create()
    # following check breaks; prob add checks to zmas in idv confs
    # kind of pain
    # ignoring the idea of storing multiple zmas right now
    # vma = automol.zmatrix.var_(zma)
    # if cnf_save_fs[0].file.vmatrix.exists():
    #     existing_vma = cnf_save_fs[0].file.vmatrix.read()
    #     assert vma == existing_vma
    # cnf_save_fs[0].file.vmatrix.write(vma)
    nsamp0 = nsamp
    inf_obj = autofile.schema.info_objects.conformer_trunk(0)
    if cnf_save_fs[0].file.info2.exists():
        inf_obj_s = cnf_save_fs[0].file.info2.read()
        nsampd = inf_obj_s.nsamp
    elif cnf_run_fs[0].file.info2.exists():
        inf_obj_r = cnf_run_fs[0].file.info2.read()
        nsampd = inf_obj_r.nsamp
    else:
        nsampd = 0

    tot_samp = nsamp - nsampd
    ioprinter.info_message(
        ' - Number of samples that have been currently run:', nsampd)
    ioprinter.info_message(' - Number of samples requested:', nsamp)

    if nsamp-nsampd > 0:
        ioprinter.info_message(
            'Running {} samples...'.format(nsamp-nsampd), newline=1)
    samp_idx = 1
    while True:
        nsamp = nsamp0 - nsampd
        # Break the while loop if enough sampls completed
        if nsamp <= 0:
            ioprinter.info_message(
                'Requested number of samples have been completed.',
                'Conformer search complete.')
            break

        # Run the conformer sampling
        if nsampd > 0:
            samp_zma, = automol.zmatrix.samples(zma, 1, tors_range_dct)
        else:
            samp_zma = zma

        ioprinter.debug_message(
            'Checking if ZMA has high repulsion...', newline=1)
        # print('zma tests:',
        #         automol.zmatrix.string(zma), automol.zmatrix.string(zma))
        bad_geom_count = 0
        while (not automol.intmol.low_repulsion_struct(zma, samp_zma) and
               bad_geom_count < 1000):
            ioprinter.warning_message('ZMA has high repulsion.', indent=1/2.)
            # print('  Bad geometry:')
            # print(automol.geom.string(automol.zmatrix.geometry(samp_zma)))
            ioprinter.warning_message(
                'Generating new sample ZMA', indent=1/2., newline=1)
            samp_zma, = automol.zmatrix.samples(zma, 1, tors_range_dct)
            bad_geom_count += 1
        ioprinter.debug_message('ZMA is fine...', indent=1/2.)

        cid = autofile.schema.generate_new_conformer_id()
        locs = [cid]

        cnf_run_fs[-1].create(locs)
        cnf_run_path = cnf_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(cnf_run_path)

        ioprinter.info_message("Run {}/{}".format(samp_idx, tot_samp))
        tors_names = list(tors_range_dct.keys())
        if two_stage and tors_names:
            ioprinter.info_message(
                'Stage one beginning, holding the coordinates constant',
                tors_names)
            es_runner.run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geom=samp_zma,
                spc_info=spc_info,
                thy_info=thy_info,
                overwrite=overwrite,
                frozen_coordinates=[tors_names],
                saddle=saddle,
                retryfail=retryfail,
                **kwargs
            )
            # print('Stage one success, reading for stage 2')
            success, ret = es_runner.read_job(
                job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if success:
                sinf_obj, _, out_str = ret
                prog = sinf_obj.prog
                samp_zma = elstruct.reader.opt_zmatrix(prog, out_str)
                ioprinter.info_message('Stage one success beginning stage two')
                # print('Stage one success beginning stage two on', samp_zma)
                es_runner.run_job(
                    job=elstruct.Job.OPTIMIZATION,
                    script_str=script_str,
                    run_fs=run_fs,
                    geom=samp_zma,
                    spc_info=spc_info,
                    thy_info=thy_info,
                    overwrite=overwrite,
                    saddle=saddle,
                    retryfail=False,
                    **kwargs
                )
        else:
            es_runner.run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geom=samp_zma,
                spc_info=spc_info,
                thy_info=thy_info,
                overwrite=overwrite,
                saddle=saddle,
                retryfail=retryfail,
                **kwargs
            )

        if cnf_save_fs[0].file.info2.exists():
            inf_obj_s = cnf_save_fs[0].file.info2.read()
            nsampd = inf_obj_s.nsamp
        elif cnf_run_fs[0].file.info2.exists():
            inf_obj_r = cnf_run_fs[0].file.info2.read()
            nsampd = inf_obj_r.nsamp
        nsampd += 1
        samp_idx += 1
        inf_obj.nsamp = nsampd
        cnf_save_fs[0].file.info2.write(inf_obj)
        cnf_run_fs[0].file.info2.write(inf_obj)


def save_conformers(cnf_run_fs, cnf_save_fs, thy_info, saddle=False,
                    rxn_class='', orig_ich=''):
    """ save the conformers that have been found so far
        # Only go through save procedure if conf not in save
        # may need to get geo, ene, etc; maybe make function
    """

    saved_locs, saved_geos, saved_enes = _saved_cnf_info(
        cnf_save_fs, thy_info)

    if not saddle:
        _check_old_inchi(orig_ich, saved_geos, saved_locs, cnf_save_fs)

    if not cnf_run_fs[0].exists():
        ioprinter.info_message(" - No conformers in run filesys to save.")
    else:
        ioprinter.info_message(" - Found conformers in run filesys to save.")
        for locs in cnf_run_fs[-1].existing():
            cnf_run_path = cnf_run_fs[-1].path(locs)
            run_fs = autofile.fs.run(cnf_run_path)
            ioprinter.reading("conformer run", cnf_run_path, newline=1)

            # Read the electronic structure optimization job
            success, ret = es_runner.read_job(
                job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)

            # Assess the geometry and save it if so
            if success:
                inf_obj, _, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                zma = elstruct.reader.opt_zmatrix(prog, out_str)
                # zma = automol.geom.zmatrix(geo)

                # Assess if geometry is properly connected
                if _geo_connected(geo, saddle):
                    if saddle:
                        ts_viable = _ts_geo_viable(
                            zma, cnf_save_fs, rxn_class, thy_info)
                        if not ts_viable:
                            continue
                    elif not _inchi_are_same(orig_ich, geo):
                        continue
                    # Assess viability of transition state conformer

                    # Determine uniqueness of conformer, save if needed
                    if _geo_unique(geo, ene, saved_geos, saved_enes, saddle):
                        # iso check breaks because of zma location
                        # if _is_proper_isomer(cnf_save_fs, zma):
                        sym_id = _sym_unique(
                            geo, ene, saved_geos, saved_enes)
                        if sym_id is None:
                            _save_unique_conformer(
                                ret, thy_info, cnf_save_fs,
                                locs, saddle=saddle)
                            saved_geos.append(geo)
                            saved_enes.append(ene)
                            saved_locs.append(locs)
                        else:
                            sym_locs = saved_locs[sym_id]
                            _save_sym_indistinct_conformer(
                                geo, cnf_save_fs, locs, sym_locs)

        # Update the conformer trajectory file
        ioprinter.obj('vspace')
        filesys.mincnf.traj_sort(cnf_save_fs, thy_info)


def _geo_connected(geo, saddle):
    """ Assess if geometry is connected. Right now only works for
        minima
    """

    # Determine connectivity (only for minima)
    if not saddle:
        gra = automol.geom.graph(geo)
        conns = automol.graph.connected_components(gra)
        lconns = len(conns)
    else:
        lconns = 1

    # Check connectivity
    if lconns == 1:
        connected = True
    else:
        ioprinter.bad_conformer('disconnected')
        connected = False

    return connected


def _geo_unique(geo, ene, seen_geos, seen_enes, saddle):
    """ Assess if a geometry is unique to saved geos
    """

    if not automol.util.value_similar_to(ene, seen_enes):
        # unique = geomprep.is_unique_tors_dist_mat_energy(
        #     geo, ene, seen_geos, seen_enes, saddle)
        unique, _ = automol.geom.compare(
                geo, seen_geos, check_dct={'dist': None, 'tors': None})
    else:
        unique = True
    if not unique:
        ioprinter.bad_conformer('not unique')

    return unique


def _inchi_are_same(orig_ich, geo):
    """ Assess if a geometry has the same connectivity to
     saved geos evaluated in temrs of inchi
    """
    same = False
    ich = automol.geom.inchi(geo)
    # if not automol.inchi.has_stereo(orig_ich):
    # if not automol.inchi.is_complete(orig_ich):
    # orig_ich = automol.inchi.add_stereo(orig_ich)
    # print(automol.inchi.has_stereo(orig_ich), 'orig_ich test', orig_ich)
    # print(automol.inchi.add_stereo(orig_ich))
    assert automol.inchi.is_complete(orig_ich), (
        'the inchi {} orig_ich is not complete'.format(orig_ich))
    if ich == orig_ich:
        same = True
    if not same:
        ioprinter.warning_message(
            " - new inchi {} not the same as old {}".format(ich, orig_ich))

    return same


def _check_old_inchi(orig_ich, seen_geos, saved_locs, cnf_save_fs):
    """
    This assumes you already have bad geos in your save
    """
    for i, geoi in enumerate(seen_geos):
        if not orig_ich == automol.geom.inchi(geoi):
            smi = automol.geom.smiles(geoi)
            ioprinter.error_message(
                'inchi do not match for {} at {}'.format(
                    smi, cnf_save_fs[-1].path(saved_locs[i])))


def _sym_unique(geo, ene, saved_geos, saved_enes, ethresh=1.0e-5):
    """ Check if a conformer is symmetrically distinct from the
        existing conformers in the filesystem
    """

    unique, sym_idx = automol.geom.compare(
        geo, seen_geos, check_dct={'coulomb': None})
    # sym_idx = None
    # for idx, (geoi, enei) in enumerate(zip(saved_geos, saved_enes)):
    #     if abs(enei - ene) < ethresh:
    #         unique = geomprep.is_unique_coulomb_energy(
    #             geo, ene, [geoi], [enei])
    #         if not unique:
    #             sym_idx = idx

    if sym_idx is not None:
        ioprinter.warning_message(' - Structure is not symmetrically unique.')

    return sym_idx


def _is_proper_isomer(cnf_save_fs, zma):
    """ Check if geom is the same isomer as those in the filesys
    """
    vma = automol.zmatrix.var_(zma)
    if cnf_save_fs[0].file.vmatrix.exists():
        exist_vma = cnf_save_fs[0].file.vmatrix.read()
        if vma != exist_vma:
            ioprinter.warning_message(
                " - Isomer is not the same as starting isomer. Skipping...")
            proper_isomer = False
        else:
            proper_isomer = True
    else:
        proper_isomer = False

    return proper_isomer


def _ts_geo_viable(zma, cnf_save_fs, rxn_class, mod_thy_info, zma_locs=(0,)):
    """ Perform a series of checks to assess the viability
        of a transition state geometry prior to saving
    """

    # Initialize viable
    viable = True

    # Obtain the min-ene zma and bond keys
    min_cnf_locs, cnf_save_path = filesys.mincnf.min_energy_conformer_locators(
        cnf_save_fs, mod_thy_info)
    zma_save_fs = fs.zmatrix(cnf_save_path)
    ref_zma = zma_save_fs[-1].file.zmatrix.read(zma_locs)

    # Read the form and broken keys from the min conf
    # frm_bnd_keys, brk_bnd_keys = tsprep.rxn_bnd_keys(
    frm_bnd_keys, brk_bnd_keys = tsprep.all_rxn_bnd_keys(
        cnf_save_fs, min_cnf_locs, zma_locs=zma_locs)

    # Use the idxs to set the forming and breaking bond names
    #    if frm_bnd_keys:
    #      frm_name = automol.zmatrix.bond_key_from_idxs(
    #           zma, frm_bnd_keys)
    #        ts_bnd1, ts_bnd2 = min(frm_bnd_keys), max(frm_bnd_keys)
    #    else:
    #      frm_name = ''
    #        ts_bnd1, ts_bnd2 = None, None

    # if brk_bnd_keys:
    #  brk_name = automol.zmatrix.bond_key_from_idxs(
    #      zma, brk_bnd_keys)
    # else:
    #  brk_name = ''
    # print('frm_name', frm_name)
    # print('brk_name', brk_name)

    # Calculate the distance of bond being formed
    # cnf_dct = automol.zmatrix.values(zma)
    # ref_dct = automol.zmatrix.values(ref_zma)
    cnf_geo = automol.zmatrix.geometry(zma)
    ref_geo = automol.zmatrix.geometry(ref_zma)

    cnf_dist_lst = []
    ref_dist_lst = []
    bnd_key_lst = []
    cnf_ang_lst = []
    ref_ang_lst = []
    for frm_bnd_key in frm_bnd_keys:
        frm_idx1, frm_idx2 = list(frm_bnd_key)
        cnf_dist = automol.geom.distance(cnf_geo, frm_idx1, frm_idx2)
        ref_dist = automol.geom.distance(ref_geo, frm_idx1, frm_idx2)
        cnf_dist_lst.append(cnf_dist)
        ref_dist_lst.append(ref_dist)
        bnd_key_lst.append(frm_bnd_key)

    for brk_bnd_key in brk_bnd_keys:
        brk_idx1, brk_idx2 = list(brk_bnd_key)
        cnf_dist = automol.geom.distance(cnf_geo, brk_idx1, brk_idx2)
        ref_dist = automol.geom.distance(ref_geo, brk_idx1, brk_idx2)
        cnf_dist_lst.append(cnf_dist)
        ref_dist_lst.append(ref_dist)
        bnd_key_lst.append(brk_bnd_key)

    for frm_bnd_key in frm_bnd_keys:
        for brk_bnd_key in brk_bnd_keys:
            for frm_idx in frm_bnd_key:
                for brk_idx in brk_bnd_key:
                    if frm_idx == brk_idx:
                        idx2 = frm_idx
                        idx1 = list(frm_bnd_key - frozenset({idx2}))[0]
                        idx3 = list(brk_bnd_key - frozenset({idx2}))[0]
                        cnf_ang = automol.geom.central_angle(
                            cnf_geo, idx1, idx2, idx3)
                        ref_ang = automol.geom.central_angle(
                            ref_geo, idx1, idx2, idx3)
                        cnf_ang_lst.append(cnf_ang)
                        ref_ang_lst.append(ref_ang)

    #      cnf_dist = cnf_dct.get(frm_name, None)
    #      ref_dist = ref_dct.get(frm_name, None)
    #  if cnf_dist is None:
    #      cnf_dist = cnf_dct.get(brk_name, None)
    #  if ref_dist is None:
    #      ref_dist = ref_dct.get(brk_name, None)
    ioprinter.debug_message('bnd_key_list', bnd_key_lst)
    ioprinter.debug_message('conf_dist', cnf_dist_lst)
    ioprinter.debug_message('ref_dist', ref_dist_lst)
    ioprinter.debug_message('conf_angle', cnf_ang_lst)
    ioprinter.debug_message('ref_angle', ref_ang_lst)

    #  # Calculate the central angle of reacting moiety of zma
    #  cnf_angle = geomprep.calc_rxn_angle(
    #      zma, frm_bnd_keys, brk_bnd_keys, rxn_class)
    #  ref_angle = geomprep.calc_rxn_angle(
    #      ref_zma, frm_bnd_keys, brk_bnd_keys, rxn_class)
    #  print('conf_angle', cnf_angle)
    #  print('ref_angle', ref_angle)

    # Set the maximum allowed displacement for a TS conformer
    max_disp = 0.6
    # would be better to check for bond forming length in bond scission with ring forming
    if 'addition' in rxn_class:
        max_disp = 0.8
    if 'abstraction' in rxn_class:
        # this was 1.4 - SJK reduced it to work for some OH abstractions
        max_disp = 1.0

    # Check forming bond angle similar to ini config
    if 'elimination' not in rxn_class:
        for ref_angle, cnf_angle in zip(ref_ang_lst, cnf_ang_lst):
            if abs(cnf_angle - ref_angle) > .44:
                ioprinter.diverged_ts('angle', ref_angle, cnf_angle)
                viable = False

    symbols = automol.geom.symbols(cnf_geo)
    lst_info = zip(ref_dist_lst, cnf_dist_lst, bnd_key_lst)
    for ref_dist, cnf_dist, bnd_key in lst_info:
        if 'add' in rxn_class or 'abst' in rxn_class:
            bnd_key1, bnd_key2 = min(list(bnd_key)), max(list(bnd_key))
            symb1 = symbols[bnd_key1]
            symb2 = symbols[bnd_key2]

            if bnd_key in frm_bnd_keys:
                # Check if radical atom is closer to some atom
                # other than the bonding atom
                cls = geomprep.is_atom_closest_to_bond_atom(
                    zma, bnd_key2, cnf_dist)
                if not cls:
                    ioprinter.diverged_ts('distance', ref_dist, cnf_dist)
                    ioprinter.info_message(
                        ' - Radical atom now has a new nearest neighbor')
                    viable = False
                # check forming bond distance
                if abs(cnf_dist - ref_dist) > max_disp:
                    ioprinter.diverged_ts('distance', ref_dist, cnf_dist)
                    viable = False

            # Check distance relative to equi. bond
            if (symb1, symb2) in bnd.LEN_DCT:
                equi_bnd = bnd.LEN_DCT[(symb1, symb2)]
            elif (symb2, symb1) in bnd.LEN_DCT:
                equi_bnd = bnd.LEN_DCT[(symb2, symb1)]
            else:
                equi_bnd = 0.0
            displace_from_equi = cnf_dist - equi_bnd
            dchk1 = abs(cnf_dist - ref_dist) > 0.1
            dchk2 = displace_from_equi < 0.2
            if dchk1 and dchk2:
                ioprinter.bad_equil_ts(cnf_dist, equi_bnd)
                viable = False
        else:
            # check forming/breaking bond distance
            # if abs(cnf_dist - ref_dist) > 0.4:
            # max disp of 0.4 causes problems for bond scission with ring forming 
            # not sure if setting it to 0.3 will cause problems for other cases
            if abs(cnf_dist - ref_dist) > 0.3:
                ioprinter.diverged_ts('distance', ref_dist, cnf_dist)
                viable = False

    return viable


def unique_fs_confs(cnf_save_fs, cnf_save_locs_lst,
                    ini_cnf_save_fs, ini_cnf_save_locs_lst):
    """ Assess which structures from the cnf_save_fs currently exist
        within the ini_cnf_save_fs. Generate a lst of unique structures
        in the ini_cnf_save_fs.
    """

    uni_ini_cnf_save_locs = []
    for ini_locs in ini_cnf_save_locs_lst:

        # Initialize variable to see if initial struct is found
        found = False

        # Loop over structs in cnf_save, see if they match the current struct
        _, inigeo = filesys.inf.cnf_fs_zma_geo(
            ini_cnf_save_fs, ini_locs)
        ini_cnf_save_path = ini_cnf_save_fs[-1].path(ini_locs)
        ioprinter.checking('structure', ini_cnf_save_path)
        for locs in cnf_save_locs_lst:
            _, geo = filesys.inf.cnf_fs_zma_geo(cnf_save_fs, locs)
            if automol.geom.almost_equal_dist_matrix(inigeo, geo, thresh=.15):
                cnf_save_path = cnf_save_fs[-1].path(locs)
                ioprinter.info_message('- Similar structure found at {}'.format(cnf_save_path))
                found = True
                break

        # If no match was found, add to unique locs lst
        if not found:
            uni_ini_cnf_save_locs.append(ini_locs)

    return uni_ini_cnf_save_locs


def _save_unique_conformer(ret, thy_info, cnf_save_fs, locs,
                           saddle=False, zma_locs=(0,)):
    """ Save the conformer in the filesystem
    """

    # Set the path to the conformer save filesystem
    cnf_save_path = cnf_save_fs[-1].path(locs)

    # Unpack the ret object and obtain the prog and method
    inf_obj, inp_str, out_str = ret
    prog = inf_obj.prog
    method = inf_obj.method

    # Read the energy and geom from the output
    ene = elstruct.reader.energy(prog, method, out_str)
    geo = elstruct.reader.opt_geometry(prog, out_str)
    zma = elstruct.reader.opt_zmatrix(prog, out_str)

    # Read the tra and graph
    if saddle:
        ts_min_cnf_locs, ts_min_path = filesys.mincnf.min_energy_conformer_locators(
            cnf_save_fs, thy_info)
        ts_min_zma_fs = fs.zmatrix(ts_min_path)
        ioprinter.debug_message('ts_min_path test:', ts_min_path)
        tra = ts_min_zma_fs[-1].file.transformation.read(zma_locs)
        ioprinter.debug_message('zma_locs test:', zma_locs)
        rct_gra = ts_min_zma_fs[-1].file.reactant_graph.read(zma_locs)

    # Build the conformer filesystem and save the structural info
    ioprinter.save_conformer(cnf_save_path)
    cnf_save_fs[-1].create(locs)
    cnf_save_fs[-1].file.geometry_info.write(inf_obj, locs)
    cnf_save_fs[-1].file.geometry_input.write(inp_str, locs)
    cnf_save_fs[-1].file.energy.write(ene, locs)
    cnf_save_fs[-1].file.geometry.write(geo, locs)

    # Build the zma filesystem and save the z-matrix
    zma_save_fs = fs.zmatrix(cnf_save_path)
    zma_save_fs[-1].create(zma_locs)
    zma_save_fs[-1].file.geometry_info.write(inf_obj, zma_locs)
    zma_save_fs[-1].file.geometry_input.write(inp_str, zma_locs)
    zma_save_fs[-1].file.zmatrix.write(zma, zma_locs)

    # Save the tra and gra for a saddle
    if saddle:
        zma_save_fs[-1].file.transformation.write(tra, zma_locs)
        zma_save_fs[-1].file.reactant_graph.write(rct_gra, zma_locs)

    # Saving the energy to a SP filesystem
    sp_save_fs = autofile.fs.single_point(cnf_save_path)
    sp_save_fs[-1].create(thy_info[1:4])
    ioprinter.save_conformer_energy(sp_save_fs[-1].root.path())
    sp_save_fs[-1].file.input.write(inp_str, thy_info[1:4])
    sp_save_fs[-1].file.info.write(inf_obj, thy_info[1:4])
    sp_save_fs[-1].file.energy.write(ene, thy_info[1:4])


def _save_sym_indistinct_conformer(geo, cnf_save_fs,
                                   cnf_tosave_locs, cnf_saved_locs):
    """ Save a structure into the SYM directory of a conformer
    """

    # Set the path to the previously saved conformer under which
    # we will save the new conformer that shares a structure
    cnf_save_path = cnf_save_fs[-1].path(cnf_saved_locs)

    # Build the sym file sys
    sym_save_fs = fs.symmetry(cnf_save_path)
    sym_save_path = cnf_save_fs[-1].path(cnf_saved_locs)
    ioprinter.save_symmetry(sym_save_path)
    sym_save_fs[-1].create(cnf_tosave_locs)
    sym_save_fs[-1].file.geometry.write(geo, cnf_tosave_locs)
    # sym_save_fs[-1].file.energy.write(ene, cnf_tosave_locs)
    # sym_save_fs[-1].file.zmatrix.write(zma, cnf_tosave_locs)


def _saved_cnf_info(cnf_save_fs, mod_thy_info):
    """ get the locs, geos and enes for saved conformers
    """

    saved_locs = list(cnf_save_fs[-1].existing())
    saved_geos = [cnf_save_fs[-1].file.geometry.read(locs)
                  for locs in saved_locs]
    saved_enes = []
    for locs in saved_locs:
        path = cnf_save_fs[-1].path(locs)
        sp_save_fs = autofile.fs.single_point(path)
        saved_enes.append(sp_save_fs[-1].file.energy.read(
            mod_thy_info[1:4]))

    return saved_locs, saved_geos, saved_enes
