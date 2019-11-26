""" drivers
"""
import os
import numpy
from scipy.interpolate import interp1d
from qcelemental import periodictable as ptab
import projrot_io
import automol
import elstruct
import autofile
import moldr
import mess_io
from submission import substr
from datalibs import phycon


def species_block(
        spc, spc_dct_i, spc_info, spc_model, pf_levels, projrot_script_str,
        elec_levels=[[0., 1]], sym_factor=1., save_prefix='spc_save_path'):
    """ prepare the species input for messpf
    """

    har_level, tors_level, vpt2_level, sym_level = pf_levels
    tors_model, vib_model, sym_model = spc_model

    # prepare the four sets of file systems
    orb_restr = moldr.util.orbital_restriction(
        spc_info, har_level)
    har_levelp = har_level[0:3]
    har_levelp.append(orb_restr)
    thy_save_fs = autofile.fs.theory(save_prefix)

    # thy_save_fs.leaf.create(har_levelp[1:4])
    har_save_path = thy_save_fs.leaf.path(har_levelp[1:4])
    saddle = False
    dist_names = []
    tors_names = []
    if 'ts_' in spc:
        har_save_fs = autofile.fs.ts(har_save_path)
        har_save_fs.trunk.create()
        har_save_path = har_save_fs.trunk.path()
        saddle = True
        tors_names = spc_dct_i['tors_names']
        if 'migration' in spc_dct_i['class'] or 'elimination' in spc_dct_i['class']:
            dist_names.append(spc_dct_i['dist_info'][0])
            dist_names.append(spc_dct_i['dist_info'][3])
    har_cnf_save_fs = autofile.fs.conformer(har_save_path)
    har_min_cnf_locs = moldr.util.min_energy_conformer_locators(har_cnf_save_fs)
    if sym_level:
        orb_restr = moldr.util.orbital_restriction(
            spc_info, sym_level)
        sym_levelp = sym_level[0:3]
        sym_levelp.append(orb_restr)

        sym_save_path = thy_save_fs.leaf.path(sym_levelp[1:4])
        if 'ts_' in spc:
            sym_save_fs = autofile.fs.ts(sym_save_path)
            sym_save_fs.trunk.create()
            sym_save_path = sym_save_fs.trunk.path()

        sym_cnf_save_fs = autofile.fs.conformer(sym_save_path)
        sym_min_cnf_locs = moldr.util.min_energy_conformer_locators(sym_cnf_save_fs)

    # Set boolean to account for a radical radical reaction (not supported by vtst)
    rad_rad_ts = False
    if 'ts_' in spc:
        if spc_dct_i['rad_rad']:
            rad_rad_ts = True

    if tors_level and not rad_rad_ts:
        orb_restr = moldr.util.orbital_restriction(
            spc_info, tors_level)
        tors_levelp = tors_level[0:3]
        tors_levelp.append(orb_restr)

        tors_save_path = thy_save_fs.leaf.path(tors_levelp[1:4])
        if 'ts_' in spc:
            tors_save_fs = autofile.fs.ts(tors_save_path)
            tors_save_fs.trunk.create()
            tors_save_path = tors_save_fs.trunk.path()
        tors_cnf_save_fs = autofile.fs.conformer(tors_save_path)
        tors_min_cnf_locs = moldr.util.min_energy_conformer_locators(tors_cnf_save_fs)
        if tors_min_cnf_locs:
            tors_cnf_save_path = tors_cnf_save_fs.leaf.path(tors_min_cnf_locs)

    if vpt2_level:
        orb_restr = moldr.util.orbital_restriction(
            spc_info, vpt2_level)
        vpt2_levelp = vpt2_level[0:3]
        vpt2_levelp.append(orb_restr)

        anh_save_path = thy_save_fs.leaf.path(vpt2_levelp[1:4])
        if 'ts_' in spc:
            anh_save_fs = autofile.fs.ts(anh_save_path)
            anh_save_fs.trunk.create()
            anh_save_path = anh_save_fs.trunk.path()

        anh_cnf_save_fs = autofile.fs.conformer(anh_save_path)
        anh_min_cnf_locs = moldr.util.min_energy_conformer_locators(anh_cnf_save_fs)
        # anh_cnf_save_path = anh_cnf_save_fs.leaf.path(anh_min_cnf_locs)

    # atom case - do as first step in each of other cases
    # pure harmonic case
    spc_str = ''
    elec_levels = [[0., spc_info[2]]]
    if 'elec_levs' in spc_dct_i:
        elec_levels = spc_dct_i['elec_levs']

    sym_factor = 1.
    form_coords = []
    if saddle:
        frm_bnd_key = spc_dct_i['frm_bnd_key']
        brk_bnd_key = spc_dct_i['brk_bnd_key']
    else:
        frm_bnd_key = []
        brk_bnd_key = []
    if 'sym' in spc_dct_i:
        sym_factor = spc_dct_i['sym']
        print('sym_factor from spc_dct_i:', sym_factor)
    else:
        if sym_model == 'SAMPLING':
            if not sym_min_cnf_locs:
                print('ERROR: Reference geometry is missing for symmetry for species {}'.format(spc_info[0]))
                return '', 0.
            sym_geo = sym_cnf_save_fs.leaf.file.geometry.read(sym_min_cnf_locs)
            sym_ene = sym_cnf_save_fs.leaf.file.energy.read(sym_min_cnf_locs)
            if dist_names:
                zma = tors_cnf_save_fs.leaf.file.zmatrix.read(tors_min_cnf_locs)
                form_coords = list(automol.zmatrix.bond_idxs(zma, dist_names[0]))
                form_coords.extend(list(dist_names[1]))
            sym_factor = moldr.conformer.symmetry_factor(
                sym_geo, sym_ene, sym_cnf_save_fs, saddle, frm_bnd_key, brk_bnd_key, form_coords, tors_names)
            print('sym_factor from conformer sampling:', sym_factor)
        if sym_model == '1DHR':
            print('Warning: the 1DHR based symmetry number has not yet been set up')
            sym_factor = 1

    imag_freq = 0.

    if (vib_model == 'HARM' and tors_model == 'RIGID') or rad_rad_ts:
        if har_min_cnf_locs is not None:
            har_geo = har_cnf_save_fs.leaf.file.geometry.read(har_min_cnf_locs)
            min_ene = har_cnf_save_fs.leaf.file.energy.read(har_min_cnf_locs)
            if automol.geom.is_atom(har_geo):
                print('This is an atom')
                mass = ptab.to_mass(har_geo[0][0])
                spc_str = mess_io.writer.atom(
                    mass, elec_levels)
            else:
                hess = har_cnf_save_fs.leaf.file.hessian.read(har_min_cnf_locs)
                freqs = elstruct.util.harmonic_frequencies(har_geo, hess, project=False)
                mode_start = 6
                if 'ts_' in spc:
                    mode_start = mode_start + 1
                    imag_freq = freqs[0]
                if automol.geom.is_linear(har_geo):
                    mode_start = mode_start - 1
                freqs = freqs[mode_start:]

                hind_rot_str = ""

                core = mess_io.writer.core_rigidrotor(har_geo, sym_factor)
                spc_str = mess_io.writer.molecule(
                    core, freqs, elec_levels,
                    hind_rot=hind_rot_str,
                    )
        else:
            print('ERROR: Reference geometry is missing for harmonic frequencies for species {}'.format(spc_info[0]))
            spc_str = ''
            imag_freq = 0.
            return '', 0.

    elif vib_model == 'HARM' and tors_model == '1DHR':
        if har_min_cnf_locs is not None:
            har_geo = har_cnf_save_fs.leaf.file.geometry.read(har_min_cnf_locs)
            min_ene = har_cnf_save_fs.leaf.file.energy.read(har_min_cnf_locs)
            if automol.geom.is_atom(har_geo):
                # print('This is an atom')
                mass = ptab.to_mass(har_geo[0][0])
                spc_str = mess_io.writer.atom(
                    mass, elec_levels)
            else:
                hess = har_cnf_save_fs.leaf.file.hessian.read(har_min_cnf_locs)
                freqs = elstruct.util.harmonic_frequencies(har_geo, hess, project=False)
                hind_rot_str = ""
                proj_rotors_str = ""
                # print('for species:', spc)
                #print('tors_min_cnf_locs test:', tors_min_cnf_locs)

                if tors_min_cnf_locs is not None:
                    if tors_cnf_save_fs.trunk.file.info.exists():
                        inf_obj_s = tors_cnf_save_fs.trunk.file.info.read()
                        tors_ranges = inf_obj_s.tors_ranges
                        tors_ranges = autofile.info.dict_(tors_ranges)
                        tors_names = list(tors_ranges.keys())
                    else:
                        print('No inf obj to identify torsional angles')
                        tors_names = []
                    zma = tors_cnf_save_fs.leaf.file.zmatrix.read(tors_min_cnf_locs)

                    tors_geo = tors_cnf_save_fs.leaf.file.geometry.read(tors_min_cnf_locs)
                    gra = automol.zmatrix.graph(zma, remove_stereo=True)
                    coo_dct = automol.zmatrix.coordinates(zma, multi=False)

                    # prepare axis, group, and projection info
                    scn_save_fs = autofile.fs.scan(tors_cnf_save_path)
                    ts_bnd = None
                    if saddle:
                        dist_name = spc_dct_i['dist_info'][0]
                        tors_names = spc_dct_i['tors_names']
                        ts_bnd = automol.zmatrix.bond_idxs(zma, dist_name)
                        ts_bnd = frozenset(ts_bnd)
                    pot = []
                    if 'hind_inc' in spc_dct_i:
                        scan_increment = spc_dct_i['hind_inc']
                    else:
                        scan_increment = 30. * phycon.DEG2RAD
                    val_dct = automol.zmatrix.values(zma)
                    tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                        zma, tors_names, scan_increment,
                        frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key)
                    tors_grids = [
                        numpy.linspace(*linspace) + val_dct[name]
                        for name, linspace in zip(tors_names, tors_linspaces)]
                    tors_sym_nums = list(automol.zmatrix.torsional_symmetry_numbers(
                        zma, tors_names, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key))
                    idx = 0
                    for tors_name, tors_grid, sym_num in zip(tors_names, tors_grids, tors_sym_nums):
                        locs_lst = []
                        enes = []
                        for grid_val in tors_grid:
                            locs_lst.append([[tors_name], [grid_val]])
                        for locs in locs_lst:
                            if scn_save_fs.leaf.exists(locs):
                                enes.append(scn_save_fs.leaf.file.energy.read(locs))
                            else:
                                enes.append(10.)
                                print('ERROR: missing grid value for torsional potential of {}'.format(spc_info[0]))

                        enes = numpy.subtract(enes, min_ene)
                        pot = list(enes*phycon.EH2KCAL)

                        # Build a potential list from only successful calculations
                        # print('pot in species block:', enes, pot)
                        pot = _hrpot_spline_fitter(pot)

                        axis = coo_dct[tors_name][1:3]

                        atm_key = axis[1]
                        if ts_bnd:
                            for atm in axis:
                                if atm in ts_bnd:
                                    atm_key = atm
                                    break
                        group = list(
                            automol.graph.branch_atom_keys(gra, atm_key, axis, saddle=saddle, ts_bnd=ts_bnd) -
                            set(axis))
                        if not group:
                            for atm in axis:
                                if atm != atm_key:
                                    atm_key = atm
                            group = list(
                                automol.graph.branch_atom_keys(gra, atm_key, axis, saddle=saddle, ts_bnd=ts_bnd) -
                                set(axis))
                        #print('sym_num before:', sym_num)
                        if saddle:
                            # check to see if fragment group was neglected
                            n_atm = automol.zmatrix.count(zma)
                            if 'addition' in spc_dct_i['class'] or 'abstraction' in spc_dct_i['class']:
                                group2 = []
                                ts_bnd1 = min(ts_bnd)
                                ts_bnd2 = max(ts_bnd)
                                for idx in range(ts_bnd2, n_atm):
                                    group2.append(idx)
                                #print('group2 test:', group2)
                                #print('axis test:', axis)
                                #print('ts_bnds test:', ts_bnd2, ts_bnd1)

                                if ts_bnd1 in group:
                                    for atm in group2:
                                        if atm not in group:
                                            group.append(atm)
                            # check to see if symmetry of XH3 rotor was missed
                            if sym_num == 1:
                                group2 = []
                                for idx in range(n_atm):
                                    if idx not in group and idx not in axis:
                                        group2.append(idx)
                                all_H = True
                                symbols = automol.zmatrix.symbols(zma)
                                #print('symbols test:', symbols)
                                #print('second group2:', group2)
                                H_count = 0
                                for idx in group2:
                                    if symbols[idx] != 'H' and symbols[idx] != 'X':
                                        all_H = False
                                        break
                                    else:
                                        if symbols[idx] == 'H':
                                            H_count += 1
                                if all_H and H_count == 3:
                                    sym_num = 3
                                    lpot = int(len(pot)/3)
                                    potp = []
                                    potp[0:lpot] = pot[0:lpot]
                                    pot = potp
                                #print('all_h test=:', all_H, H_count)
                                
                        #print('sym_num after:', sym_num)

                        group = list(numpy.add(group, 1))
                        axis = list(numpy.add(axis, 1))
                        #print('axis test:', axis)
                        #print('atm_key:', atm_key)
                        #print('group:', group)
                        #for idx, atm in enumerate(axis):
                        #    if atm == atm_key+1:
                        #        if idx != 0:
                        #            axis.reverse()
                        #            print('axis reversed', axis)
                        if (atm_key+1) != axis[1]:
                            axis.reverse()
                            #print('axis reversed:', axis)
                        #if atm_key != axis(0):
                            #axis.reverse()

                        #check for dummy transformations
                        atom_symbols = automol.zmatrix.symbols(zma)
                        dummy_idx = []
                        for atm_idx, atm in enumerate(atom_symbols):
                            if atm == 'X':
                                dummy_idx.append(atm_idx)
                        remdummy = numpy.zeros(len(zma[0]))
                        for dummy in dummy_idx:
                            for idx, _ in enumerate(remdummy):
                                if dummy < idx:
                                   remdummy[idx] += 1
                        hind_rot_str += mess_io.writer.rotor_hindered(
                            group, axis, sym_num, pot, remdummy=remdummy)
                        #print('projrot 1 test:')
                        proj_rotors_str += projrot_io.writer.rotors(
                            axis, group, remdummy=remdummy)
                        sym_factor /= sym_num
                        idx += 1

                    # create a messpf input file and run messpf to get tors_freqs and tors_zpes
                    if saddle and tors_names is not None:
                        dummy_freqs = [1000.]
                        dummy_zpe = 0.0
                        core = mess_io.writer.core_rigidrotor(tors_geo, sym_factor)
                        spc_str = mess_io.writer.molecule(
                            core, dummy_freqs, elec_levels,
                            hind_rot=hind_rot_str,
                            )
                        temp_step = 100.
                        ntemps = 5
                        zpe_str = '{0:<8.2f}\n'.format(dummy_zpe)
                        zpe_str = ' ZeroEnergy[kcal/mol] ' + zpe_str
                        zpe_str += 'End\n'
                        global_pf_str = mess_io.writer.global_pf(
                            [], temp_step, ntemps, rel_temp_inc=0.001,
                            atom_dist_min=0.6)
                        spc_head_str = 'Species ' + ' Tmp'
                        pf_inp_str = '\n'.join(
                            [global_pf_str, spc_head_str,
                             spc_str, zpe_str])

                        bld_locs = ['PF', 0]
                        bld_save_fs = autofile.fs.build(tors_save_path)
                        bld_save_fs.leaf.create(bld_locs)
                        pf_path = bld_save_fs.leaf.path(bld_locs)

                        # run messpf
                        with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
                            pf_file.write(pf_inp_str)
                        pf_script_str = ("#!/usr/bin/env bash\n"
                                         "export OMP_NUM_THREADS=10\n"
                                         "messpf pf.inp pf.out >> stdout.log &> stderr.log")

                        moldr.util.run_script(pf_script_str, pf_path)

                        with open(os.path.join(pf_path, 'pf.log'), 'r') as mess_file:
                            output_string = mess_file.read()

                        # Read the freqs and zpes
                        tors_freqs = mess_io.reader.tors.freqs(output_string)
                        tors_zpes = mess_io.reader.tors.zpves(output_string)
                    else:
                        tors_freqs = []
                        tors_zpes = []

                    tors_zpe = 0.0
                    for (tors_freq, tors_1dhr_zpe) in zip(tors_freqs, tors_zpes):
                        tors_zpe += tors_1dhr_zpe

                    # run one version of ProjRot to get the projected frequencies for that version
                    # Write the string for the ProjRot input
                    coord_proj = 'cartesian'
                    grad = ''
                    #print('projrot 2 test:')
                    projrot_inp_str = projrot_io.writer.rpht_input(
                        tors_geo, grad, hess, rotors_str=proj_rotors_str,
                        coord_proj=coord_proj)

                    bld_locs = ['PROJROT', 0]
                    bld_save_fs = autofile.fs.build(tors_save_path)
                    bld_save_fs.leaf.create(bld_locs)
                    path = bld_save_fs.leaf.path(bld_locs)
                    print('Build Path for Partition Functions in species block')
                    print(path)
                    proj_file_path = os.path.join(path, 'RPHt_input_data.dat')
                    with open(proj_file_path, 'w') as proj_file:
                        proj_file.write(projrot_inp_str)

                    moldr.util.run_script(projrot_script_str, path)

                    freqs = []
                    zpe_har_no_tors = 0.
                    har_zpe = 0.
                    if pot:
                        rthrproj_freqs, _ = projrot_io.reader.rpht_output(
                            path+'/hrproj_freq.dat')
                        freqs = rthrproj_freqs
                        zpe_har_no_tors = sum(freqs)*phycon.WAVEN2KCAL/2.
                    rtproj_freqs, imag_freq = projrot_io.reader.rpht_output(
                        path+'/RTproj_freq.dat')
                    har_zpe = sum(rtproj_freqs)*phycon.WAVEN2KCAL/2.
                    if not freqs:
                        freqs = rtproj_freqs
                    if 'ts_' in spc:
                        if imag_freq:
                            imag_freq = imag_freq[0]
                        else:
                            imag_freq = freqs[-1]
                            freqs = freqs[:-1]

                    # now run the other version of ProjRot
                    projrot_script_str2 = ("#!/usr/bin/env bash\n"
                    "RPHt.exe >& /dev/null")
                    moldr.util.run_script(projrot_script_str2, path)
                    zpe_har_no_tors_2 = 0.0
                    freqs_2 = []
                    if pot:
                        rthrproj_freqs_2, _ = projrot_io.reader.rpht_output(
                            path+'/hrproj_freq.dat')
                        freqs_2 = rthrproj_freqs_2
                        zpe_har_no_tors_2 = sum(freqs_2)*phycon.WAVEN2KCAL/2.
                    rtproj_freqs, imag_freq_2 = projrot_io.reader.rpht_output(
                        path+'/RTproj_freq.dat')
                    har_zpe = sum(rtproj_freqs)*phycon.WAVEN2KCAL/2.
                    if not freqs_2:
                        freqs_2 = rtproj_freqs
                    if 'ts_' in spc:
                        if imag_freq_2:
                            imag_freq_2 = imag_freq_2[0]
                        else:
                            imag_freq_2 = freqs_2[-1]
                            freqs_2 = freqs_2[:-1]

                    har_tors_zpe = har_zpe - zpe_har_no_tors
                    har_tors_zpe_2 = har_zpe - zpe_har_no_tors_2
                    del_tors_zpe = har_tors_zpe - tors_zpe
                    del_tors_zpe_2 = har_tors_zpe_2 - tors_zpe
                    # print('tors_zpe test:', del_tors_zpe, del_tors_zpe_2)
                    if del_tors_zpe <= del_tors_zpe_2:
                        zpe = zpe_har_no_tors + tors_zpe
                    else:
                        zpe = zpe_har_no_tors_2 + tors_zpe
                        freqs = freqs_2
                        imag_freq = imag_freq_2
                    # shutil.rmtree(path)
                    # print('freqs test in species block', freqs, imag_freq)
                    # print('zpe test in species_block:',zpe_har_no_tors, zpe_har_no_tors_2, 
                    #        tors_zpe, har_zpe, zpe)
                    core = mess_io.writer.core_rigidrotor(tors_geo, sym_factor)
                    spc_str = mess_io.writer.molecule(
                        core, freqs, elec_levels,
                        hind_rot=hind_rot_str
                        )
        else:
            print('ERROR: Reference geometry is missing for harmonic frequencies for species {}'
                  .format(spc_info[0]))
            spc_str = ''
            imag_freq = 0.
            return '', 0.

    elif vib_model == 'HARM' and tors_model == 'MDHR':
        print('HARM and MDHR combination is not yet implemented')

    elif vib_model == 'HARM' and tors_model == 'TAU':
        print('HARM and TAU combination is not yet implemented')
        moldr.driver.tau_pf_write(
            name=name,
            save_prefix=thy_save_path,
            run_grad=run_grad_pf,
            run_hess=run_hess_pf,
        )
    elif vib_model == 'VPT2' and tors_model == 'RIGID':
        if anh_min_cnf_locs is not None:
            anh_geo = anh_cnf_save_fs.leaf.file.geometry.read(anh_min_cnf_locs)
            min_ene = anh_cnf_save_fs.leaf.file.energy.read(anh_min_cnf_locs)
            if automol.geom.is_atom(anh_geo):
                # print('This is an atom')
                mass = ptab.to_mass(anh_geo[0][0])
                spc_str = mess_io.writer.atom(
                    mass, elec_levels)
            else:
                hess = anh_cnf_save_fs.leaf.file.hessian.read(anh_min_cnf_locs)
                freqs = elstruct.util.harmonic_frequencies(anh_geo, hess, project=True)
                mode_start = 6
                if 'ts_' in spc:
                    mode_start = mode_start + 1
                    imag_freq = freqs[0]
                if automol.geom.is_linear(anh_geo):
                    mode_start = mode_start - 1
                freqs = freqs[mode_start:]
                zpe = sum(freqs)*phycon.WAVEN2KCAL/2.
                hind_rot_str = ""

                core = mess_io.writer.core_rigidrotor(anh_geo, sym_factor)
                spc_str = mess_io.writer.molecule(
                    core, freqs, elec_levels,
                    hind_rot=hind_rot_str,
                    )
        else:
            print('ERROR: Reference geometry is missing for anharmonic analysis for species {}'.format(spc_info[0]))
            spc_str = ''
            imag_freq = 0.
            return '', 0.
        print('VPT2 and RIGID combination is not yet properly implemented')

    elif vib_model == 'VPT2' and tors_model == '1DHR':
        print('VPT2 and 1DHR combination is not yet implemented')

    elif vib_model == 'VPT2' and tors_model == 'TAU':
        print('VPT2 and TAU combination is not yet implemented')
        moldr.driver.tau_pf_write(
            name=name,
            save_prefix=thy_save_path,
            run_grad=run_grad_pf,
            run_hess=run_hess_pf,
        )

    return spc_str, imag_freq


def vtst_with_no_saddle_block(
        ts_dct, ts_label, reac_label, prod_label, spc_ene, rct_zpe, projrot_script_str,
        multi_info, elec_levels=[[0., 1]], sym_factor=1.
        ):
    """ prepare the mess input string for a variational TS that does not have
    a saddle point. Do it by calling the species block for each grid point
    in the scan file system
    """

    ts_info = ['', ts_dct['chg'], ts_dct['mul']]
    print('ts_dct test:', ts_dct['mul'])
    print('multi info test:', multi_info)
    orb_restr = moldr.util.orbital_restriction(ts_info, multi_info)
    multi_level = multi_info[0:3]
    multi_level.append(orb_restr)

    rxn_run_path = ts_dct['rxn_fs'][2]
    thy_run_fs = autofile.fs.theory(rxn_run_path)
    thy_run_fs.leaf.create(multi_level[1:4])
    thy_run_path = thy_run_fs.leaf.path(multi_level[1:4])

    rxn_save_path = ts_dct['rxn_fs'][3]
    thy_save_fs = autofile.fs.theory(rxn_save_path)
    thy_save_fs.leaf.create(multi_level[1:4])
    thy_save_path = thy_save_fs.leaf.path(multi_level[1:4])

    scn_run_fs = autofile.fs.scan(thy_run_path)
    scn_save_fs = autofile.fs.scan(thy_save_path)

    # read the scan save file system to get the energies, zero-point energies,
    # geometries, hessians,
    # ultimately

    sym_factor = 1.
    irc_pt_strs = []
    proj_rotors_str = ''
    coord_proj = 'cartesian'
    grad = ''
    pot = []

    elec_levels=[[0., ts_dct['mul']]]
    grid = ts_dct['grid']
    grid1 = grid[0]
    grid2 = grid[1]
    grid = numpy.append(grid[0], grid[1])
    dist_name = ts_dct['dist_info'][0]

    inf_locs = [[dist_name], [1000.]]
    inf_sep_ene = scn_save_fs.leaf.file.energy.read(inf_locs)
    # print('inf sep ene in vtst with no saddle:', inf_sep_ene)
    # print('grid test in vtst with no saddle:', grid)

    grid[::-1].sort()
    for idx, grid_val in enumerate(grid):
        # print('idx, grid_val test:', idx, grid_val)
        locs = [[dist_name], [grid_val]]
        # print('scn save fs:', scn_save_fs.leaf.path(locs))

        # get geometry
        if not scn_save_fs.leaf.file.geometry.exists(locs):
            continue
        else:
            geom = scn_save_fs.leaf.file.geometry.read(locs)

        # get energy
        if not scn_save_fs.leaf.file.energy.exists(locs):
            continue
        else:
            ene = scn_save_fs.leaf.file.energy.read(locs)

        # get gradient
        #if not scn_save_fs.leaf.file.gradient.exists(locs):
            #continue
        #else:
            #grad = scn_save_fs.leaf.file.gradient.read(locs)
            #print('grad in vtst: \n', grad)

        # get hessian
        if not scn_save_fs.leaf.file.hessian.exists(locs):
            continue
        else:
            hess = scn_save_fs.leaf.file.hessian.read(locs)

            projrot_inp_str = projrot_io.writer.rpht_input(
                geom, grad, hess, rotors_str=proj_rotors_str,
                coord_proj=coord_proj)

            scn_save_path= scn_save_fs.leaf.path(locs)
            bld_locs = ['PROJROT', 0]
            bld_save_fs = autofile.fs.build(scn_save_path)
            bld_save_fs.leaf.create(bld_locs)
            path = bld_save_fs.leaf.path(bld_locs)
            print('Build Path for Partition Functions')
            print(path)
            proj_file_path = os.path.join(path, 'RPHt_input_data.dat')
            with open(proj_file_path, 'w') as proj_file:
                proj_file.write(projrot_inp_str)

            moldr.util.run_script(projrot_script_str, path)

            freqs = []
            if len(pot) > 0:
                rthrproj_freqs, _ = projrot_io.reader.rpht_output(
                    path+'/hrproj_freq.dat')
                freqs = rthrproj_freqs
            rtproj_freqs, imag_freq = projrot_io.reader.rpht_output(
                path+'/RTproj_freq.dat')
            if not freqs:
                freqs = rtproj_freqs
                if not imag_freq:
                    freqs = freqs[:-1]

            freqs_test_0 = elstruct.util.harmonic_frequencies(geom, hess, project=False)
            mode_start = 7
            if automol.geom.is_linear(geom):
                mode_start = mode_start - 1
            freqs_test = freqs_test_0[mode_start:]
            print('projrot freqs in vtst:', freqs)
            #print('all unprojected freqs in vtst:', freqs_test_0)
            #print('unprojected freqs in vtst:', freqs_test)

        zpe = sum(freqs)*phycon.WAVEN2KCAL/2.
        # for now use the zpe calculated at the first grid point as an approximation to
        # the zpe at infinite sepration
        if idx == 0:
            rct_zpe = zpe

        erel = (ene - inf_sep_ene)*phycon.EH2KCAL
        erel_zpe_corr = erel + zpe - rct_zpe
        #eref = erel
        eref_abs = erel_zpe_corr + spc_ene
        # print('vtst inf test:', ene, inf_sep_ene)
        # print('vtst ene test:', eref_abs, erel_zpe_corr, erel, zpe, spc_ene, rct_zpe)

        # Iniialize the header of the string
        irc_pt_str = '!----------------------------------------------- \n'
        irc_pt_str += '! IRC Point {0}\n'.format(str(idx+1))

        # Write the molecule section for each irc point
        core = mess_io.writer.mol_data.core_rigidrotor(geom, sym_factor, interp_emax='')
        irc_pt_str += mess_io.writer.species.molecule(core, freqs, elec_levels,
             hind_rot='', xmat=None, rovib_coups='', rot_dists='')

        # Append the zero energy for the molecule
        irc_pt_str += '    ZeroEnergy[kcal/mol]      {0:<8.2f}\n'.format(eref_abs)
        if grid_val != grid[-1]:
            irc_pt_str += 'End \n'

        # Append string to list
        irc_pt_strs.append(irc_pt_str)

    # Write the MESS string for the variational sections
    variational_str = mess_io.writer.rxnchan.ts_variational(
        ts_label, reac_label, prod_label, irc_pt_strs)
    # print('variational_str test:', variational_str)

    return variational_str


def vtst_saddle_block(scn_save_fs, geoms, frequencies, energies):
    """ prepare the mess input string for a variational TS where there is a 
    saddle point on the MEP. In this case, there is limited torsional information.
    """

    # read the scan save file system to get the energies, zero-point energies, symmetry numbers, 
    # geometries, hessians, torsional potentials for each point on the MEP

    # Determine the the number of points along the irc
    nirc = 21

    # Loop over all the points of the irc and build MESS strings
    irc_pt_strings = []
    for i in range(nirc):

        # Iniialize the header of the string
        irc_pt_string = '!-----------------------------------------------'
        irc_pt_string += '! IRC Point {0}\n'.format(str(i+1))

        # Write the molecule section for each irc point
        core = mess_io.writer.mol_data.core_rigidrotor(geom1, sym_factor, interp_emax='')
        irc_pt_str += mess_io.writer.species.molecule(core, freqs, elec_levels,
             hind_rot='', xmat=None, rovib_coups='', rot_dists='')

        # Append the zero point energy for the molecule
        irc_pt_str += '    ZeroEnergy[kcal/mol]      {0:<8.2f}'.format(zero_energy)

        # Append string to list
        irc_pt_strings.append(irc_pt_string)

    # Write the MESS string for the variational sections
    varational_str = mess_io.writer.rxnchan.ts_variational(
        ts_label, reac_label, prod_label, irc_pt_strings)

    return variational_str


def pst_block(
        spc_dct_i, spc_dct_j, spc_model, pf_levels, projrot_script_str,
        spc_save_fs, elec_levels=[[0., 1]], sym_factor=1.,
        pst_params=[1.0, 6]
        ):
    """ prepare a Phase Space Theory species block
    """

    har_level, tors_level, vpt2_level, sym_level = pf_levels
    tors_model, vib_model, sym_model = spc_model

    # prepare the four sets of file systems
    spc_info_i = (spc_dct_i['ich'], spc_dct_i['chg'], spc_dct_i['mul'])
    spc_info_j = (spc_dct_j['ich'], spc_dct_j['chg'], spc_dct_j['mul'])
    spc_save_fs.leaf.create(spc_info_i)
    spc_save_fs.leaf.create(spc_info_j)
    save_path_i = spc_save_fs.leaf.path(spc_info_i)
    save_path_j = spc_save_fs.leaf.path(spc_info_j)

    orb_restr = moldr.util.orbital_restriction(
        spc_info_i, har_level)
    har_levelp_i = har_level[0:3]
    har_levelp_i.append(orb_restr)
    thy_save_fs_i = autofile.fs.theory(save_path_i)
    orb_restr = moldr.util.orbital_restriction(
        spc_info_j, har_level)
    har_levelp_j = har_level[0:3]
    har_levelp_j.append(orb_restr)
    thy_save_fs_j = autofile.fs.theory(save_path_j)

    har_save_path_i = thy_save_fs_i.leaf.path(har_levelp_i[1:4])
    har_save_path_j = thy_save_fs_j.leaf.path(har_levelp_j[1:4])
    har_cnf_save_fs_i = autofile.fs.conformer(har_save_path_i)
    har_cnf_save_fs_j = autofile.fs.conformer(har_save_path_j)
    har_min_cnf_locs_i = moldr.util.min_energy_conformer_locators(har_cnf_save_fs_i)
    har_min_cnf_locs_j = moldr.util.min_energy_conformer_locators(har_cnf_save_fs_j)

    if sym_level:
        orb_restr = moldr.util.orbital_restriction(
            spc_info_i, sym_level)
        sym_levelp_i = sym_level[0:3]
        sym_levelp_i.append(orb_restr)
        sym_save_path_i = thy_save_fs_i.leaf.path(sym_levelp_i[1:4])
        sym_cnf_save_fs_i = autofile.fs.conformer(sym_save_path_i)
        sym_min_cnf_locs_i = moldr.util.min_energy_conformer_locators(sym_cnf_save_fs_i)
        orb_restr = moldr.util.orbital_restriction(
            spc_info_j, sym_level)
        sym_levelp_j = sym_level[0:3]
        sym_levelp_j.append(orb_restr)
        sym_save_path_j = thy_save_fs_j.leaf.path(sym_levelp_j[1:4])
        sym_cnf_save_fs_j = autofile.fs.conformer(sym_save_path_j)
        sym_min_cnf_locs_j = moldr.util.min_energy_conformer_locators(sym_cnf_save_fs_j)

    if tors_level:
        orb_restr = moldr.util.orbital_restriction(
            spc_info_i, tors_level)
        tors_levelp_i = tors_level[0:3]
        tors_levelp_i.append(orb_restr)
        tors_save_path_i = thy_save_fs_i.leaf.path(tors_levelp_i[1:4])
        tors_cnf_save_fs_i = autofile.fs.conformer(tors_save_path_i)
        tors_min_cnf_locs_i = moldr.util.min_energy_conformer_locators(tors_cnf_save_fs_i)
        tors_cnf_save_path_i = tors_cnf_save_fs_i.leaf.path(tors_min_cnf_locs_i)
        orb_restr = moldr.util.orbital_restriction(
            spc_info_j, tors_level)
        tors_levelp_j = tors_level[0:3]
        tors_levelp_j.append(orb_restr)
        tors_save_path_j = thy_save_fs_j.leaf.path(tors_levelp_j[1:4])
        tors_cnf_save_fs_j = autofile.fs.conformer(tors_save_path_j)
        tors_min_cnf_locs_j = moldr.util.min_energy_conformer_locators(tors_cnf_save_fs_j)
        tors_cnf_save_path_j = tors_cnf_save_fs_j.leaf.path(tors_min_cnf_locs_j)

    if vpt2_level:
        orb_restr = moldr.util.orbital_restriction(
            spc_info_i, vpt2_level)
        vpt2_levelp_i = vpt2_level[0:3]
        vpt2_levelp_i.append(orb_restr)
        anh_save_path_i = thy_save_fs_i.leaf.path(vpt2_levelp_i[1:4])
        anh_cnf_save_fs_i = autofile.fs.conformer(anh_save_path_i)
        anh_min_cnf_locs_i = moldr.util.min_energy_conformer_locators(anh_cnf_save_fs_i)
        orb_restr = moldr.util.orbital_restriction(
            spc_info_j, vpt2_level)
        vpt2_levelp_j = vpt2_level[0:3]
        vpt2_levelp_j.append(orb_restr)
        anh_save_path_j = thy_save_fs_j.leaf.path(vpt2_levelp_j[1:4])
        anh_cnf_save_fs_j = autofile.fs.conformer(anh_save_path_j)
        anh_min_cnf_locs_j = moldr.util.min_energy_conformer_locators(anh_cnf_save_fs_j)

    spc_str = ''
    if 'elec_levs' in spc_dct_i:
        elec_levels_i = spc_dct_i['elec_levs']
    else:
        elec_levels_i = [[0., spc_dct_i['mul']]]
    if 'elec_levs' in spc_dct_j:
        elec_levels_j = spc_dct_j['elec_levs']
    else:
        elec_levels_j = [[0., spc_dct_j['mul']]]

    # Combine the energy levels
    init_elec_levels = []
    for _, elec_level_i in enumerate(elec_levels_i):
        for _, elec_level_j in enumerate(elec_levels_j):
            init_elec_levels.append(
                [elec_level_i[0]+elec_level_j[0],
                 elec_level_i[1]*elec_level_j[1]])

    # See if any levels repeat and thus need to be added together
    elec_levels = []
    for level in init_elec_levels:
        # Put level in in final list
        if level not in elec_levels:
            elec_levels.append(level)
        # Add the level to the one in the list
        else:
            idx = elec_levels.index(level)
            elec_levels[idx][1] += level[1]

    sym_factor_i = 1.
    sym_factor_j = 1.
    if 'sym' in spc_dct_i:
        sym_factor_i = spc_dct_i['sym']
    else:
        if sym_model == 'SAMPLING':
            sym_geo_i = sym_cnf_save_fs_i.leaf.file.geometry.read(sym_min_cnf_locs_i)
            sym_ene_i = sym_cnf_save_fs_i.leaf.file.energy.read(sym_min_cnf_locs_i)
            sym_factor_i = moldr.conformer.symmetry_factor(sym_geo_i, sym_ene_i, sym_cnf_save_fs_i)
        if sym_model == '1DHR':
            # Warning: the 1DHR based symmetry number has not yet been set up
            sym_factor_i = 1
    if 'sym' in spc_dct_j:
        sym_factor_j = spc_dct_j['sym']
    else:
        if sym_model == 'SAMPLING':
            sym_geo_j = sym_cnf_save_fs_j.leaf.file.geometry.read(sym_min_cnf_locs_j)
            sym_ene_j = sym_cnf_save_fs_j.leaf.file.energy.read(sym_min_cnf_locs_j)
            sym_factor_j = moldr.conformer.symmetry_factor(sym_geo_j, sym_ene_j, sym_cnf_save_fs_j)
        if sym_model == '1DHR':
            # Warning: the 1DHR based symmetry number has not yet been set up
            sym_factor_j = 1
    sym_factor = sym_factor_i * sym_factor_j

    if vib_model == 'HARM' and tors_model == 'RIGID':
        if har_min_cnf_locs_i is not None:
            har_geo_i = har_cnf_save_fs_i.leaf.file.geometry.read(har_min_cnf_locs_i)
            if har_min_cnf_locs_j is not None:
                har_geo_j = har_cnf_save_fs_j.leaf.file.geometry.read(har_min_cnf_locs_j)
                freqs = []
                freqs_i = []
                freqs_j = []
                if not automol.geom.is_atom(har_geo_i):
                    hess_i = har_cnf_save_fs_i.leaf.file.hessian.read(har_min_cnf_locs_i)
                    freqs_i = elstruct.util.harmonic_frequencies(har_geo_i, hess_i, project=False)
                    mode_start = 6
                    if automol.geom.is_linear(har_geo_i):
                        mode_start = mode_start - 1
                    freqs = freqs_i[mode_start:]
                if not automol.geom.is_atom(har_geo_j):
                    hess_j = har_cnf_save_fs_j.leaf.file.hessian.read(har_min_cnf_locs_j)
                    freqs_j = elstruct.util.harmonic_frequencies(
                        har_geo_j, hess_j, project=False)
                    mode_start = 6
                    if automol.geom.is_linear(har_geo_j):
                        mode_start = mode_start - 1
                    freqs += freqs_j[mode_start:]
                hind_rot_str = ""
                form_i = automol.geom.formula(har_geo_i)
                form_j = automol.geom.formula(har_geo_j)
                form = automol.formula.join(form_i, form_j)
                stoich = ''
                for key, val in form.items():
                    stoich += key + str(val)
                core = mess_io.writer.core_phasespace(
                    har_geo_i, har_geo_j, sym_factor, stoich,
                    pot_prefactor=pst_params[0], pot_power_exp=pst_params[1])
                spc_str = mess_io.writer.molecule(
                    core, freqs, elec_levels,
                    hind_rot=hind_rot_str,
                    )
        else:
            spc_str = ''

    if vib_model == 'HARM' and tors_model == '1DHR':
        if har_min_cnf_locs_i is not None:
            har_geo_i = har_cnf_save_fs_i.leaf.file.geometry.read(har_min_cnf_locs_i)
            min_ene_i = har_cnf_save_fs_i.leaf.file.energy.read(har_min_cnf_locs_i)
            if har_min_cnf_locs_j is not None:
                har_geo_j = har_cnf_save_fs_j.leaf.file.geometry.read(har_min_cnf_locs_j)
                min_ene_j = har_cnf_save_fs_j.leaf.file.energy.read(har_min_cnf_locs_j)

                freqs_i = []
                freqs_j = []
                is_atom_i = automol.geom.is_atom(har_geo_i)
                is_atom_j = automol.geom.is_atom(har_geo_j)
                if not is_atom_i:
                    hess_i = har_cnf_save_fs_i.leaf.file.hessian.read(har_min_cnf_locs_i)
                    freqs_i = elstruct.util.harmonic_frequencies(har_geo_i, hess_i, project=False)
                    mode_start = 6
                    if automol.geom.is_linear(har_geo_i):
                        mode_start = mode_start - 1
                    freqs_i = freqs_i[mode_start:]
                if not is_atom_j:
                    hess_j = har_cnf_save_fs_j.leaf.file.hessian.read(har_min_cnf_locs_j)
                    freqs_j = elstruct.util.harmonic_frequencies(har_geo_j, hess_j, project=False)
                    mode_start = 6
                    if automol.geom.is_linear(har_geo_j):
                        mode_start = mode_start - 1
                    freqs_j = freqs_j[mode_start:]

                proj_rotors_str = ""
                hind_rot_str = ""

                if tors_min_cnf_locs_i is not None and not is_atom_i:
                    if har_cnf_save_fs_i.trunk.file.info.exists():
                        inf_obj_s = har_cnf_save_fs_i.trunk.file.info.read()
                        tors_ranges = inf_obj_s.tors_ranges
                        tors_ranges = autofile.info.dict_(tors_ranges)
                        tors_names = list(tors_ranges.keys())
                    else:
                        print('No inf obj to identify torsional angles')
                        tors_names = []
                    zma = tors_cnf_save_fs_i.leaf.file.zmatrix.read(tors_min_cnf_locs_i)

                    tors_geo = tors_cnf_save_fs_i.leaf.file.geometry.read(tors_min_cnf_locs_i)
                    gra = automol.zmatrix.graph(zma, remove_stereo=True)
                    coo_dct = automol.zmatrix.coordinates(zma, multi=False)

                    # prepare axis, group, and projection info
                    scn_save_fs = autofile.fs.scan(tors_cnf_save_path_i)
                    pot = []
                    if 'hind_inc' in spc_dct_i:
                        scan_increment = spc_dct_i['hind_inc']
                    else:
                        scan_increment = 30. * phycon.DEG2RAD
                    val_dct = automol.zmatrix.values(zma)
                    tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                        zma, tors_names, scan_increment)
                    tors_grids = [
                        numpy.linspace(*linspace) + val_dct[name]
                        for name, linspace in zip(tors_names, tors_linspaces)]
                    tors_sym_nums = list(automol.zmatrix.torsional_symmetry_numbers(
                        zma, tors_names))
                    for tors_name, tors_grid, sym_num in zip(tors_names, tors_grids, tors_sym_nums):
                        locs_lst = []
                        enes = []
                        for grid_val in tors_grid:
                            locs_lst.append([[tors_name], [grid_val]])
                        for locs in locs_lst:
                            if scn_save_fs.leaf.exists(locs):
                                enes.append(scn_save_fs.leaf.file.energy.read(locs))
                            else:
                                enes.append(10.)
                                print('ERROR: missing grid value for torsional potential of {}'
                                      .format(spc_info_i[0]))
                        enes = numpy.subtract(enes, min_ene_i)
                        pot = list(enes*phycon.EH2KCAL)

                        # Build a potential list from only successful calculations
                        pot = _hrpot_spline_fitter(pot)

                        axis = coo_dct[tors_name][1:3]

                        atm_key = axis[1]
                        group = list(
                            automol.graph.branch_atom_keys(gra, atm_key, axis) - set(axis))
                        if not group:
                            for atm in axis:
                                if atm != atm_key:
                                    atm_key = atm
                            group = list(
                                automol.graph.branch_atom_keys(gra, atm_key, axis) - set(axis))

                        group = list(numpy.add(group, 1))
                        axis = list(numpy.add(axis, 1))
                        #print('axis test:', axis)
                        #print('atm_key:', atm_key)
                        #print('group:', group)
                        #for idx, atm in enumerate(axis):
                        #    if atm == atm_key+1:
                        #        if idx != 0:
                        #            axis.reverse()
                        #            print('axis reversed', axis)
                        if (atm_key+1) != axis[1]:
                            axis.reverse()
                            #print('axis reversed:', axis)
                        #if atm_key != axis(0):
                            #axis.reverse()

                        #check for dummy transformations
                        atom_symbols = automol.zmatrix.symbols(zma)
                        dummy_idx = []
                        for atm_idx, atm in enumerate(atom_symbols):
                            if atm == 'X':
                                dummy_idx.append(atm_idx)
                        remdummy = numpy.zeros(len(zma[0]))
                        for dummy in dummy_idx:
                            for idx, _ in enumerate(remdummy):
                                if dummy < idx:
                                    remdummy[idx] += 1
                        hind_rot_str += mess_io.writer.rotor_hindered(
                            group, axis, sym_num, pot, remdummy=remdummy, geom=har_geo_i)
                        #print('projrot 3 test:')
                        proj_rotors_str += projrot_io.writer.rotors(
                            axis, group, remdummy=remdummy)
                        sym_factor /= sym_num

                    # Write the string for the ProjRot input
                    coord_proj = 'cartesian'
                    grad = ''
                    #print('projrot 4 test:')
                    projrot_inp_str = projrot_io.writer.rpht_input(
                        tors_geo, grad, hess_i, rotors_str=proj_rotors_str,
                        coord_proj=coord_proj)

                    bld_locs = ['PROJROT', 0]
                    bld_save_fs = autofile.fs.build(tors_save_path_i)
                    bld_save_fs.leaf.create(bld_locs)
                    path = bld_save_fs.leaf.path(bld_locs)
                    print('Build Path for Partition Functions')
                    print(path)
                    proj_file_path = os.path.join(path, 'RPHt_input_data.dat')
                    with open(proj_file_path, 'w') as proj_file:
                        proj_file.write(projrot_inp_str)

                    moldr.util.run_script(projrot_script_str, path)

                    freqs_i = []
                    if pot:
                        rthrproj_freqs, _ = projrot_io.reader.rpht_output(
                            path+'/hrproj_freq.dat')
                        freqs_i = rthrproj_freqs
                    if not freqs_i:
                        rtproj_freqs, _ = projrot_io.reader.rpht_output(
                            path+'/RTproj_freq.dat')
                        freqs_i = rtproj_freqs

                proj_rotors_str = ""
                if tors_min_cnf_locs_j is not None and not is_atom_j:
                    if har_cnf_save_fs_j.trunk.file.info.exists():
                        inf_obj_s = har_cnf_save_fs_j.trunk.file.info.read()
                        tors_ranges = inf_obj_s.tors_ranges
                        tors_ranges = autofile.info.dict_(tors_ranges)
                        tors_names = list(tors_ranges.keys())
                    else:
                        print('No inf obj to identify torsional angles')
                        tors_names = []
                    zma = tors_cnf_save_fs_j.leaf.file.zmatrix.read(tors_min_cnf_locs_j)

                    tors_geo = tors_cnf_save_fs_j.leaf.file.geometry.read(tors_min_cnf_locs_j)
                    gra = automol.zmatrix.graph(zma, remove_stereo=True)
                    coo_dct = automol.zmatrix.coordinates(zma, multi=False)

                    # prepare axis, group, and projection info
                    scn_save_fs = autofile.fs.scan(tors_cnf_save_path_j)
                    pot = []
                    if 'hind_inc' in spc_dct_j:
                        scan_increment = spc_dct_j['hind_inc']
                    else:
                        scan_increment = 30. * phycon.DEG2RAD
                    val_dct = automol.zmatrix.values(zma)
                    tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                        zma, tors_names, scan_increment)
                    tors_grids = [
                        numpy.linspace(*linspace) + val_dct[name]
                        for name, linspace in zip(tors_names, tors_linspaces)]
                    tors_sym_nums = list(automol.zmatrix.torsional_symmetry_numbers(
                        zma, tors_names))
                    for tors_name, tors_grid, sym_num in zip(tors_names, tors_grids, tors_sym_nums):
                        locs_lst = []
                        enes = []
                        for grid_val in tors_grid:
                            locs_lst.append([[tors_name], [grid_val]])
                        for locs in locs_lst:
                            if scn_save_fs.leaf.exists(locs):
                                enes.append(scn_save_fs.leaf.file.energy.read(locs))
                            else:
                                enes.append(10.)
                                print('ERROR: missing grid value for torsional potential of {}'
                                      .format(spc_info_j[0]))
                        enes = numpy.subtract(enes, min_ene_j)
                        pot = list(enes*phycon.EH2KCAL)

                        # Build a potential list from only successful calculations
                        pot = _hrpot_spline_fitter(pot)

                        axis = coo_dct[tors_name][1:3]

                        atm_key = axis[1]
                        group = list(
                            automol.graph.branch_atom_keys(gra, atm_key, axis) - set(axis))
                        if not group:
                            for atm in axis:
                                if atm != atm_key:
                                    atm_key = atm
                            group = list(
                                automol.graph.branch_atom_keys(gra, atm_key, axis) - set(axis))

                        group = list(numpy.add(group, 1))
                        axis = list(numpy.add(axis, 1))
                        #print('axis test:', axis)
                        #print('atm_key:', atm_key)
                        #print('group:', group)
                        if (atm_key+1) != axis[1]:
                            axis.reverse()
                            #print('axis reversed:', axis)

                        #check for dummy transformations
                        atom_symbols = automol.zmatrix.symbols(zma)
                        dummy_idx = []
                        for atm_idx, atm in enumerate(atom_symbols):
                            if atm == 'X':
                                dummy_idx.append(atm_idx)
                        remdummy = numpy.zeros(len(zma[0]))
                        for dummy in dummy_idx:
                            for idx, _ in enumerate(remdummy):
                                if dummy < idx:
                                    remdummy[idx] += 1
                        hind_rot_str += mess_io.writer.rotor_hindered(
                            group, axis, sym_num, pot, remdummy, geom=har_geo_j)
                        #print('projrot 5 test:')
                        proj_rotors_str += projrot_io.writer.rotors(
                            axis, group, remdummy=remdummy)
                        sym_factor /= sym_num

                    # Write the string for the ProjRot input
                    coord_proj = 'cartesian'
                    grad = ''
                    #print('projrot 6 test:')
                    projrot_inp_str = projrot_io.writer.rpht_input(
                        tors_geo, grad, hess_j, rotors_str=proj_rotors_str,
                        coord_proj=coord_proj)

                    bld_locs = ['PROJROT', 0]
                    bld_save_fs = autofile.fs.build(tors_save_path_j)
                    bld_save_fs.leaf.create(bld_locs)
                    path = bld_save_fs.leaf.path(bld_locs)
                    print('Build Path for Partition Functions')
                    print(path)
                    proj_file_path = os.path.join(path, 'RPHt_input_data.dat')
                    with open(proj_file_path, 'w') as proj_file:
                        proj_file.write(projrot_inp_str)

                    moldr.util.run_script(projrot_script_str, path)

                    freqs_j = []
                    if pot:
                        rthrproj_freqs, _ = projrot_io.reader.rpht_output(
                            path+'/hrproj_freq.dat')
                        freqs_j = rthrproj_freqs
                    if not freqs_j:
                        rtproj_freqs, imag_freq = projrot_io.reader.rpht_output(
                            path+'/RTproj_freq.dat')
                        freqs_j = rtproj_freqs

                freqs = list(freqs_i) + list(freqs_j)

                form_i = automol.geom.formula(har_geo_i)
                form_j = automol.geom.formula(har_geo_j)
                form = automol.formula.join(form_i, form_j)
                stoich = ''
                for key, val in form.items():
                    stoich += key + str(val)
                core = mess_io.writer.core_phasespace(
                    har_geo_i, har_geo_j, sym_factor, stoich,
                    pot_prefactor=pst_params[0], pot_power_exp=pst_params[1])
                spc_str = mess_io.writer.molecule(
                    core, freqs, elec_levels,
                    hind_rot=hind_rot_str,
                    )

        else:
            spc_str = ''

    return spc_str


def fake_species_block(
        spc_dct_i, spc_dct_j, spc_info_i, spc_info_j, spc_model, pf_levels, projrot_script_str,
        elec_levels=[[0., 1]], sym_factor=1.,
        save_prefix_i='spc_save_path', save_prefix_j='spc_save_path'):
    """ prepare a fake species block corresponding to the van der Waals well between two fragments
    """
    har_level, tors_level, vpt2_level, sym_level = pf_levels
    tors_model, vib_model, sym_model = spc_model

    # prepare the four sets of file systems
    orb_restr = moldr.util.orbital_restriction(
        spc_info_i, har_level)
    har_levelp_i = har_level[0:3]
    har_levelp_i.append(orb_restr)
    thy_save_fs_i = autofile.fs.theory(save_prefix_i)
    orb_restr = moldr.util.orbital_restriction(
        spc_info_j, har_level)
    har_levelp_j = har_level[0:3]
    har_levelp_j.append(orb_restr)
    thy_save_fs_j = autofile.fs.theory(save_prefix_j)

    # thy_save_fs.leaf.create(har_levelp[1:4])
    har_save_path_i = thy_save_fs_i.leaf.path(har_levelp_i[1:4])
    har_save_path_j = thy_save_fs_j.leaf.path(har_levelp_j[1:4])
    har_cnf_save_fs_i = autofile.fs.conformer(har_save_path_i)
    har_cnf_save_fs_j = autofile.fs.conformer(har_save_path_j)
    har_min_cnf_locs_i = moldr.util.min_energy_conformer_locators(har_cnf_save_fs_i)
    har_min_cnf_locs_j = moldr.util.min_energy_conformer_locators(har_cnf_save_fs_j)

    if sym_level:
        orb_restr = moldr.util.orbital_restriction(
            spc_info_i, sym_level)
        sym_levelp_i = sym_level[0:3]
        sym_levelp_i.append(orb_restr)
        sym_save_path_i = thy_save_fs_i.leaf.path(sym_levelp_i[1:4])
        sym_cnf_save_fs_i = autofile.fs.conformer(sym_save_path_i)
        sym_min_cnf_locs_i = moldr.util.min_energy_conformer_locators(sym_cnf_save_fs_i)

        orb_restr = moldr.util.orbital_restriction(
            spc_info_j, sym_level)
        sym_levelp_j = sym_level[0:3]
        sym_levelp_j.append(orb_restr)
        sym_save_path_j = thy_save_fs_j.leaf.path(sym_levelp_j[1:4])
        sym_cnf_save_fs_j = autofile.fs.conformer(sym_save_path_j)
        sym_min_cnf_locs_j = moldr.util.min_energy_conformer_locators(sym_cnf_save_fs_j)

    tors_names = []
    if tors_level:
        orb_restr = moldr.util.orbital_restriction(
            spc_info_i, tors_level)
        tors_levelp_i = tors_level[0:3]
        tors_levelp_i.append(orb_restr)
        orb_restr = moldr.util.orbital_restriction(
            spc_info_j, tors_level)
        tors_levelp_j = tors_level[0:3]
        tors_levelp_j.append(orb_restr)

        tors_save_path_i = thy_save_fs_i.leaf.path(tors_levelp_i[1:4])
        tors_cnf_save_fs_i = autofile.fs.conformer(tors_save_path_i)
        tors_min_cnf_locs_i = moldr.util.min_energy_conformer_locators(tors_cnf_save_fs_i)
        tors_cnf_save_path_i = tors_cnf_save_fs_i.leaf.path(tors_min_cnf_locs_i)

        tors_save_path_j = thy_save_fs_j.leaf.path(tors_levelp_j[1:4])
        tors_cnf_save_fs_j = autofile.fs.conformer(tors_save_path_j)
        tors_min_cnf_locs_j = moldr.util.min_energy_conformer_locators(tors_cnf_save_fs_j)
        tors_cnf_save_path_j = tors_cnf_save_fs_j.leaf.path(tors_min_cnf_locs_j)

    spc_str = ''
    if 'elec_levs' in spc_dct_i:
        elec_levels_i = spc_dct_i['elec_levs']
    else:
        elec_levels_i = [[0., spc_dct_i['mul']]]
    if 'elec_levs' in spc_dct_j:
        elec_levels_j = spc_dct_j['elec_levs']
    else:
        elec_levels_j = [[0., spc_dct_j['mul']]]

    # Combine the energy levels
    init_elec_levels = []
    for _, elec_level_i in enumerate(elec_levels_i):
        for _, elec_level_j in enumerate(elec_levels_j):
            init_elec_levels.append(
                [elec_level_i[0]+elec_level_j[0],
                 elec_level_i[1]*elec_level_j[1]])

    # See if any levels repeat and thus need to be added together
    elec_levels = []
    for level in init_elec_levels:
        # Put level in in final list
        if level not in elec_levels:
            elec_levels.append(level)
        # Add the level to the one in the list
        else:
            idx = elec_levels.index(level)
            elec_levels[idx][1] += level[1]

    # print('elec_levels final ')
    # print(elec_levels)

    sym_factor = 1.
    sym_factor_i = 1.
    sym_factor_j = 1.
    if 'sym' in spc_dct_i:
        sym_factor_i = spc_dct_i['sym']
    else:
        if sym_model == 'SAMPLING':
            sym_geo_i = sym_cnf_save_fs_i.leaf.file.geometry.read(sym_min_cnf_locs_i)
            sym_ene_i = sym_cnf_save_fs_i.leaf.file.energy.read(sym_min_cnf_locs_i)
            sym_factor_i = moldr.conformer.symmetry_factor(sym_geo_i, sym_ene_i, sym_cnf_save_fs_i)
        if sym_model == '1DHR':
            # Warning: the 1DHR based symmetry number has not yet been set up
            sym_factor_i = 1
    if 'sym' in spc_dct_j:
        sym_factor_j = spc_dct_j['sym']
    else:
        if sym_model == 'SAMPLING':
            sym_geo_j = sym_cnf_save_fs_j.leaf.file.geometry.read(sym_min_cnf_locs_j)
            sym_ene_j = sym_cnf_save_fs_j.leaf.file.energy.read(sym_min_cnf_locs_j)
            sym_factor_j = moldr.conformer.symmetry_factor(sym_geo_j, sym_ene_j, sym_cnf_save_fs_j)
        if sym_model == '1DHR':
            # Warning: the 1DHR based symmetry number has not yet been set up
            sym_factor_j = 1
    sym_factor = sym_factor_i * sym_factor_j

    if vib_model == 'HARM' and tors_model == 'RIGID':
        if har_min_cnf_locs_i is not None:
            har_geo_i = har_cnf_save_fs_i.leaf.file.geometry.read(har_min_cnf_locs_i)
            if har_min_cnf_locs_j is not None:
                har_geo_j = har_cnf_save_fs_j.leaf.file.geometry.read(har_min_cnf_locs_j)

                freqs = [30, 50, 70, 100, 200]
                ntrans = 5
                is_atom_i = automol.geom.is_atom(har_geo_i)
                is_linear_i = automol.geom.is_linear(har_geo_i)
                is_atom_j = automol.geom.is_atom(har_geo_j)
                is_linear_j = automol.geom.is_linear(har_geo_i)
                if is_atom_i:
                    ntrans = ntrans - 3
                if is_atom_j:
                    ntrans = ntrans - 3
                if is_linear_i:
                    ntrans = ntrans - 2
                if is_linear_j:
                    ntrans = ntrans - 2
                if is_atom_i and is_atom_j:
                    ntrans = 0
                freqs = freqs[0:ntrans]
                if not is_atom_i:
                    hess_i = har_cnf_save_fs_i.leaf.file.hessian.read(har_min_cnf_locs_i)
                    freqs_i = elstruct.util.harmonic_frequencies(har_geo_i, hess_i, project=False)
                    mode_start = 6
                    if automol.geom.is_linear(har_geo_i):
                        mode_start = mode_start - 1
                    freqs += freqs_i[mode_start:]
                if not is_atom_j:
                    hess_j = har_cnf_save_fs_j.leaf.file.hessian.read(har_min_cnf_locs_j)
                    freqs_j = elstruct.util.harmonic_frequencies(har_geo_j, hess_j, project=False)
                    mode_start = 6
                    if automol.geom.is_linear(har_geo_j):
                        mode_start = mode_start - 1
                    freqs += freqs_j[mode_start:]

                max_z_i = max(atom[1][2] for atom in har_geo_i)
                min_z_j = min(atom[1][2] for atom in har_geo_j)
                har_geo = har_geo_i
                har_geo_j = automol.geom.translated(har_geo_j, [0., 0., max_z_i + 6. - min_z_j])
                har_geo += har_geo_j

                hind_rot_str = ""

                core = mess_io.writer.core_rigidrotor(har_geo, sym_factor)
                spc_str = mess_io.writer.molecule(
                    core, freqs, elec_levels,
                    hind_rot=hind_rot_str,
                    )
        else:
            spc_str = ''

    if vib_model == 'HARM' and tors_model == '1DHR':
        if har_min_cnf_locs_i is not None:
            har_geo_i = har_cnf_save_fs_i.leaf.file.geometry.read(har_min_cnf_locs_i)
            min_ene_i = har_cnf_save_fs_i.leaf.file.energy.read(har_min_cnf_locs_i)
            if har_min_cnf_locs_j is not None:
                har_geo_j = har_cnf_save_fs_j.leaf.file.geometry.read(har_min_cnf_locs_j)
                min_ene_j = har_cnf_save_fs_j.leaf.file.energy.read(har_min_cnf_locs_j)
                har_geo_js = har_geo_j

                freqs_trans = [30, 50, 70, 100, 200]
                ntrans = 5
                is_atom_i = automol.geom.is_atom(har_geo_i)
                is_linear_i = automol.geom.is_linear(har_geo_i)
                is_atom_j = automol.geom.is_atom(har_geo_j)
                is_linear_j = automol.geom.is_linear(har_geo_i)
                if is_atom_i:
                    ntrans = ntrans - 3
                if is_atom_j:
                    ntrans = ntrans - 3
                if is_linear_i:
                    ntrans = ntrans - 2
                if is_linear_j:
                    ntrans = ntrans - 2
                if is_atom_i and is_atom_j:
                    ntrans = 0
                freqs_trans = freqs_trans[0:ntrans]
                freqs_i = []
                freqs_j = []
                if not is_atom_i:
                    hess_i = har_cnf_save_fs_i.leaf.file.hessian.read(har_min_cnf_locs_i)
                    freqs_i = elstruct.util.harmonic_frequencies(har_geo_i, hess_i, project=False)
                    mode_start = 6
                    if automol.geom.is_linear(har_geo_i):
                        mode_start = mode_start - 1
                    freqs_i = freqs_i[mode_start:]
                if not is_atom_j:
                    hess_j = har_cnf_save_fs_j.leaf.file.hessian.read(har_min_cnf_locs_j)
                    freqs_j = elstruct.util.harmonic_frequencies(har_geo_j, hess_j, project=False)
                    mode_start = 6
                    if automol.geom.is_linear(har_geo_j):
                        mode_start = mode_start - 1
                    freqs_j = freqs_j[mode_start:]

                max_z_i = max(atom[1][2] for atom in har_geo_i)
                min_z_j = min(atom[1][2] for atom in har_geo_j)
                har_geo = har_geo_i
                har_geo_j = automol.geom.translated(har_geo_j, [0., 0., max_z_i + 6. - min_z_j])
                har_geo += har_geo_j

                hind_rot_str = ""
                proj_rotors_str = ""
                #print('fb tors_min_cnf_locs_i test:', tors_min_cnf_locs_i)

                if tors_min_cnf_locs_i is not None and not is_atom_i:
                    if har_cnf_save_fs_i.trunk.file.info.exists():
                        inf_obj_s = har_cnf_save_fs_i.trunk.file.info.read()
                        tors_ranges = inf_obj_s.tors_ranges
                        tors_ranges = autofile.info.dict_(tors_ranges)
                        tors_names = list(tors_ranges.keys())
                    else:
                        print('No inf obj to identify torsional angles')
                        tors_names = []
                    zma = tors_cnf_save_fs_i.leaf.file.zmatrix.read(tors_min_cnf_locs_i)

                    tors_geo = tors_cnf_save_fs_i.leaf.file.geometry.read(tors_min_cnf_locs_i)
                    gra = automol.zmatrix.graph(zma, remove_stereo=True)
                    coo_dct = automol.zmatrix.coordinates(zma, multi=False)

                    # prepare axis, group, and projection info
                    scn_save_fs = autofile.fs.scan(tors_cnf_save_path_i)
                    pot = []
                    if 'hind_inc' in spc_dct_i:
                        scan_increment = spc_dct_i['hind_inc']
                    else:
                        scan_increment = 30. * phycon.DEG2RAD
                    val_dct = automol.zmatrix.values(zma)
                    tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                        zma, tors_names, scan_increment)
                    tors_grids = [
                        numpy.linspace(*linspace) + val_dct[name]
                        for name, linspace in zip(tors_names, tors_linspaces)]
                    tors_sym_nums = list(automol.zmatrix.torsional_symmetry_numbers(
                        zma, tors_names))
                    #print('fb tors_names test:', tors_names)
                    for tors_name, tors_grid, sym_num in zip(tors_names, tors_grids, tors_sym_nums):
                        locs_lst = []
                        enes = []
                        for grid_val in tors_grid:
                            locs_lst.append([[tors_name], [grid_val]])
                        for locs in locs_lst:
                            if scn_save_fs.leaf.exists(locs):
                                enes.append(scn_save_fs.leaf.file.energy.read(locs))
                            else:
                                enes.append(10.)
                                print('ERROR: missing grid value for torsional potential of {}'
                                      .format(spc_info_i[0]))
                        enes = numpy.subtract(enes, min_ene_i)
                        pot = list(enes*phycon.EH2KCAL)

                        # Build a potential list from only successful calculations
                        pot = _hrpot_spline_fitter(pot)
                        # print('fb pot test:', pot)

                        axis = coo_dct[tors_name][1:3]

                        atm_key = axis[1]
                        group = list(
                            automol.graph.branch_atom_keys(gra, atm_key, axis) - set(axis))
                        if not group:
                            for atm in axis:
                                if atm != atm_key:
                                    atm_key = atm
                            group = list(
                                automol.graph.branch_atom_keys(gra, atm_key, axis) - set(axis))

                        group = list(numpy.add(group, 1))
                        axis = list(numpy.add(axis, 1))
                        #print('fb axis test:', axis)
                        #print('fb atm_key:', atm_key)
                        #print('fb group:', group)
                        #for idx, atm in enumerate(axis):
                        #    if atm == atm_key+1:
                        #        if idx != 0:
                        #            axis.reverse()
                        #            print('axis reversed', axis)
                        if (atm_key+1) != axis[1]:
                            axis.reverse()
                            #print('fb axis reversed:', axis)
                        #if atm_key != axis(0):
                            #axis.reverse()

                        #check for dummy transformations
                        atom_symbols = automol.zmatrix.symbols(zma)
                        dummy_idx = []
                        for atm_idx, atm in enumerate(atom_symbols):
                            if atm == 'X':
                                dummy_idx.append(atm_idx)
                        remdummy = numpy.zeros(len(zma[0]))
                        for dummy in dummy_idx:
                            for idx, _ in enumerate(remdummy):
                                if dummy < idx:
                                    remdummy[idx] += 1
                        hind_rot_str += mess_io.writer.rotor_hindered(
                            group, axis, sym_num, pot, remdummy=remdummy, geom=har_geo_i)
                        proj_rotors_str += projrot_io.writer.rotors(
                            axis, group, remdummy=remdummy)
                        #print('projrot 7 test:', proj_rotors_str)
                        sym_factor /= sym_num

                    # Write the string for the ProjRot input
                    coord_proj = 'cartesian'
                    grad = ''
                    #print('projrot 8 test:', proj_rotors_str)
                    projrot_inp_str = projrot_io.writer.rpht_input(
                        tors_geo, grad, hess_i, rotors_str=proj_rotors_str,
                        coord_proj=coord_proj)

                    bld_locs = ['PROJROT', 0]
                    bld_save_fs = autofile.fs.build(tors_save_path_i)
                    bld_save_fs.leaf.create(bld_locs)
                    path = bld_save_fs.leaf.path(bld_locs)
                    print('Build Path for Partition Functions')
                    print(path)
                    proj_file_path = os.path.join(path, 'RPHt_input_data.dat')
                    with open(proj_file_path, 'w') as proj_file:
                        proj_file.write(projrot_inp_str)

                    moldr.util.run_script(projrot_script_str, path)

                    freqs_i = []
                    if pot:
                        rthrproj_freqs, _ = projrot_io.reader.rpht_output(
                            path+'/hrproj_freq.dat')
                        freqs_i = rthrproj_freqs
                    if not freqs_i:
                        rtproj_freqs, _ = projrot_io.reader.rpht_output(
                            path+'/RTproj_freq.dat')
                        freqs_i = rtproj_freqs

                proj_rotors_str = ""
                if tors_min_cnf_locs_j is not None and not is_atom_j:
                    if har_cnf_save_fs_j.trunk.file.info.exists():
                        inf_obj_s = har_cnf_save_fs_j.trunk.file.info.read()
                        tors_ranges = inf_obj_s.tors_ranges
                        tors_ranges = autofile.info.dict_(tors_ranges)
                        tors_names = list(tors_ranges.keys())
                    else:
                        print('No inf obj to identify torsional angles')
                        tors_names = []
                    zma = tors_cnf_save_fs_j.leaf.file.zmatrix.read(tors_min_cnf_locs_j)

                    tors_geo = tors_cnf_save_fs_j.leaf.file.geometry.read(tors_min_cnf_locs_j)
                    gra = automol.zmatrix.graph(zma, remove_stereo=True)
                    coo_dct = automol.zmatrix.coordinates(zma, multi=False)

                    # prepare axis, group, and projection info
                    scn_save_fs = autofile.fs.scan(tors_cnf_save_path_j)
                    pot = []
                    if 'hind_inc' in spc_dct_j:
                        scan_increment = spc_dct_j['hind_inc']
                    else:
                        scan_increment = 30. * phycon.DEG2RAD
                    val_dct = automol.zmatrix.values(zma)
                    tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                        zma, tors_names, scan_increment)
                    tors_grids = [
                        numpy.linspace(*linspace) + val_dct[name]
                        for name, linspace in zip(tors_names, tors_linspaces)]
                    tors_sym_nums = list(automol.zmatrix.torsional_symmetry_numbers(
                        zma, tors_names))
                    for tors_name, tors_grid, sym_num in zip(tors_names, tors_grids, tors_sym_nums):
                        locs_lst = []
                        enes = []
                        for grid_val in tors_grid:
                            locs_lst.append([[tors_name], [grid_val]])
                        for locs in locs_lst:
                            if scn_save_fs.leaf.exists(locs):
                                enes.append(scn_save_fs.leaf.file.energy.read(locs))
                            else:
                                enes.append(10.)
                                print('ERROR: missing grid value for torsional potential of {}'
                                      .format(spc_info_j[0]))
                        enes = numpy.subtract(enes, min_ene_j)
                        pot = list(enes*phycon.EH2KCAL)

                        # Build a potential list from only successful calculations
                        pot = _hrpot_spline_fitter(pot)

                        axis = coo_dct[tors_name][1:3]

                        atm_key = axis[1]
                        group = list(
                            automol.graph.branch_atom_keys(gra, atm_key, axis) - set(axis))
                        if not group:
                            for atm in axis:
                                if atm != atm_key:
                                    atm_key = atm
                            group = list(
                                automol.graph.branch_atom_keys(gra, atm_key, axis) - set(axis))

                        group = list(numpy.add(group, 1))
                        axis = list(numpy.add(axis, 1))
                        #print('axis test:', axis)
                        #print('atm_key:', atm_key)
                        #print('group:', group)
                        if (atm_key+1) != axis[1]:
                            axis.reverse()
                            #print('axis reversed:', axis)

                        #check for dummy transformations
                        atom_symbols = automol.zmatrix.symbols(zma)
                        dummy_idx = []
                        for atm_idx, atm in enumerate(atom_symbols):
                            if atm == 'X':
                                dummy_idx.append(atm_idx)
                        remdummy = numpy.zeros(len(zma[0]))
                        for dummy in dummy_idx:
                            for idx, _ in enumerate(remdummy):
                                if dummy < idx:
                                    remdummy[idx] += 1
                        hind_rot_str += mess_io.writer.rotor_hindered(
                            group, axis, sym_num, pot, remdummy=remdummy, geom=har_geo_js)
                        #print('projrot 9 test:')
                        proj_rotors_str += projrot_io.writer.rotors(
                            axis, group, remdummy=remdummy)
                        sym_factor /= sym_num

                    # Write the string for the ProjRot input
                    coord_proj = 'cartesian'
                    grad = ''
                    #print('projrot 10 test:')
                    projrot_inp_str = projrot_io.writer.rpht_input(
                        tors_geo, grad, hess_j, rotors_str=proj_rotors_str,
                        coord_proj=coord_proj)

                    bld_locs = ['PROJROT', 0]
                    bld_save_fs = autofile.fs.build(tors_save_path_j)
                    bld_save_fs.leaf.create(bld_locs)
                    path = bld_save_fs.leaf.path(bld_locs)
                    print('Build Path for Partition Functions')
                    print(path)
                    proj_file_path = os.path.join(path, 'RPHt_input_data.dat')
                    with open(proj_file_path, 'w') as proj_file:
                        proj_file.write(projrot_inp_str)

                    moldr.util.run_script(projrot_script_str, path)

                    freqs_j = []
                    if pot:
                        rthrproj_freqs, _ = projrot_io.reader.rpht_output(
                            path+'/hrproj_freq.dat')
                        freqs_j = rthrproj_freqs
                    if not freqs_j:
                        rtproj_freqs, imag_freq = projrot_io.reader.rpht_output(
                            path+'/RTproj_freq.dat')
                        freqs_j = rtproj_freqs

                #print('freq_trans test:', freqs_trans)
                #print('freq_i test:', list(freqs_i))
                #print('freq_j test:', list(freqs_j))
                freqs = freqs_trans + list(freqs_i) + list(freqs_j)

                core = mess_io.writer.core_rigidrotor(har_geo, sym_factor)
                spc_str = mess_io.writer.molecule(
                    core, freqs, elec_levels,
                    hind_rot=hind_rot_str,
                    )
        else:
            spc_str = ''
    return spc_str


def get_high_level_energy(
        spc_info, thy_low_level, thy_high_level, save_prefix, saddle=False):
    """ get high level energy at low level optimized geometry
    """
    if saddle:
        spc_save_path = save_prefix
    else:
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs.leaf.create(spc_info)
        spc_save_path = spc_save_fs.leaf.path(spc_info)

    orb_restr = moldr.util.orbital_restriction(
        spc_info, thy_low_level)
    thy_low_level = thy_low_level[1:3]
    thy_low_level.append(orb_restr)

    ll_save_fs = autofile.fs.theory(spc_save_path)
    ll_save_path = ll_save_fs.leaf.path(thy_low_level)

    if saddle:
        ll_save_fs = autofile.fs.ts(ll_save_path)
        ll_save_fs.trunk.create()
        ll_save_path = ll_save_fs.trunk.path()

    cnf_save_fs = autofile.fs.conformer(ll_save_path)
    min_cnf_locs = moldr.util.min_energy_conformer_locators(
        cnf_save_fs)
    if not min_cnf_locs:
        print('ERROR: No minimum conformer geometry for this species {}'.format(spc_info[0]))
        return 0.0
    cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
    # min_cnf_geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)

    orb_restr = moldr.util.orbital_restriction(
        spc_info, thy_high_level)
    thy_high_level = thy_high_level[1:3]
    thy_high_level.append(orb_restr)

    sp_save_fs = autofile.fs.single_point(cnf_save_path)
    sp_save_fs.leaf.create(thy_high_level)

    min_ene = sp_save_fs.leaf.file.energy.read(thy_high_level)

    return min_ene


def get_zero_point_energy(
        spc, spc_dct_i, pf_levels, spc_model, pf_script_str,
        elec_levels=[[0., 1]], sym_factor=1.,
        save_prefix='spc_save_path'):
    """ compute the ZPE including torsional and anharmonic corrections
    """

    projrot_script_str = substr.PROJROT

    spc_info = (spc_dct_i['ich'], spc_dct_i['chg'], spc_dct_i['mul'])
    # prepare the sets of file systems
    har_level, tors_level, vpt2_level, _ = pf_levels
    tors_model, vib_model, _ = spc_model

    thy_save_fs = autofile.fs.theory(save_prefix)

    orb_restr = moldr.util.orbital_restriction(
        spc_info, har_level)
    har_levelp = har_level[0:3]
    har_levelp.append(orb_restr)

    har_save_path = thy_save_fs.leaf.path(har_levelp[1:4])
    saddle = False
    if 'ts_' in spc:
        har_save_fs = autofile.fs.ts(har_save_path)
        har_save_fs.trunk.create()
        har_save_path = har_save_fs.trunk.path()
        saddle = True

    # print('inside zpe saddle is:', saddle)
    har_cnf_save_fs = autofile.fs.conformer(har_save_path)
    har_min_cnf_locs = moldr.util.min_energy_conformer_locators(har_cnf_save_fs)
    
    # Set boolean to account for a radical radical reaction (not supported by vtst)
    rad_rad_ts = False
    if 'ts_' in spc:
        if spc_dct_i['rad_rad']:
            rad_rad_ts = True

    if tors_level and not rad_rad_ts:
        orb_restr = moldr.util.orbital_restriction(
            spc_info, tors_level)
        tors_levelp = tors_level[0:3]
        tors_levelp.append(orb_restr)

        # thy_save_fs.leaf.create(tors_levelp[1:4])
        tors_save_path = thy_save_fs.leaf.path(tors_levelp[1:4])
        if 'ts_' in spc:
            tors_save_fs = autofile.fs.ts(tors_save_path)
            tors_save_fs.trunk.create()
            tors_save_path = tors_save_fs.trunk.path()

        tors_cnf_save_fs = autofile.fs.conformer(tors_save_path)
        tors_min_cnf_locs = moldr.util.min_energy_conformer_locators(tors_cnf_save_fs)
        # print('tors_save_path test:', tors_save_path)
        # print('tors_min_cnf_locs test:', tors_min_cnf_locs)
        if tors_min_cnf_locs:
            tors_cnf_save_path = tors_cnf_save_fs.leaf.path(tors_min_cnf_locs)

    if vpt2_level:
        orb_restr = moldr.util.orbital_restriction(
            spc_info, vpt2_level)
        vpt2_levelp = vpt2_level[0:3]
        vpt2_levelp.append(orb_restr)

        anh_save_path = thy_save_fs.leaf.path(vpt2_levelp[1:4])
        if 'ts_' in spc:
            anh_save_fs = autofile.fs.ts(anh_save_path)
            anh_save_fs.trunk.create()
            anh_save_path = anh_save_fs.trunk.path()

        anh_cnf_save_fs = autofile.fs.conformer(anh_save_path)
        anh_min_cnf_locs = moldr.util.min_energy_conformer_locators(anh_cnf_save_fs)

    if saddle:
        frm_bnd_key = spc_dct_i['frm_bnd_key']
        brk_bnd_key = spc_dct_i['brk_bnd_key']
    else:
        frm_bnd_key = []
        brk_bnd_key = []
    har_zpe = 0.0
    is_atom = False
    # get reference harmonic
    if not har_min_cnf_locs:
        print('ERROR: No harmonic reference geometry for this species {}'.format(spc_info[0]))
        return har_zpe, is_atom
    har_geo = har_cnf_save_fs.leaf.file.geometry.read(har_min_cnf_locs)
    if automol.geom.is_atom(har_geo):
        har_zpe = 0.0
        is_atom = True

    else:
        hess = har_cnf_save_fs.leaf.file.hessian.read(har_min_cnf_locs)
        freqs = elstruct.util.harmonic_frequencies(har_geo, hess, project=False)

        mode_start = 6
        if 'ts_' in spc:
            mode_start = mode_start + 1
        if automol.geom.is_linear(har_geo):
            mode_start = mode_start - 1
        freqs = freqs[mode_start:]

        har_zpe = sum(freqs)*phycon.WAVEN2KCAL/2.

    # if 'radical radical' in spc['class']:
    #    ret = har_zpe

    print('vib_model in zpe:', vib_model)
    print('tors_model in zpe:', tors_model)
    if (vib_model == 'HARM' and tors_model == 'RIGID') or rad_rad_ts:
        ret = har_zpe

    elif vib_model == 'HARM' and tors_model == '1DHR':
        # make pf string for 1d rotor
        # run messpf
        # read 1d harmonic and torsional ZPEs
        # modify har_zpe

        zpe = har_zpe
        hind_rot_str = ""
        proj_rotors_str = ""
        # print('tors_min_cnf_locs:', tors_min_cnf_locs)
        tors_names = []
        if tors_min_cnf_locs is not None:
            if tors_cnf_save_fs.trunk.file.info.exists():
                inf_obj_s = har_cnf_save_fs.trunk.file.info.read()
                tors_ranges = inf_obj_s.tors_ranges
                tors_ranges = autofile.info.dict_(tors_ranges)
                #print(tors_ranges)
                tors_names = list(tors_ranges.keys())

            else:
                print('No inf obj to identify torsional angles')

            min_ene = tors_cnf_save_fs.leaf.file.energy.read(tors_min_cnf_locs)
            tors_geo = tors_cnf_save_fs.leaf.file.geometry.read(tors_min_cnf_locs)
            zma = tors_cnf_save_fs.leaf.file.zmatrix.read(tors_min_cnf_locs)
            gra = automol.zmatrix.graph(zma, remove_stereo=True)
            tors_zpe_cor = 0.0
            coo_dct = automol.zmatrix.coordinates(zma, multi=False)
            # prepare axis, group, info
            scn_save_fs = autofile.fs.scan(tors_cnf_save_path)
            ts_bnd = None
            if saddle:
                dist_name = spc_dct_i['dist_info'][0]
                tors_names = spc_dct_i['tors_names']
                ts_bnd =  automol.zmatrix.bond_idxs(zma, dist_name)
                ts_bnd = frozenset(ts_bnd)
            pot = []
            if 'hind_inc' in spc_dct_i:
                scan_increment = spc_dct_i['hind_inc']
            else:
                scan_increment = 30. * phycon.DEG2RAD
            val_dct = automol.zmatrix.values(zma)
            tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                zma, tors_names, scan_increment, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key)
            tors_grids = [numpy.linspace(*linspace) + val_dct[name]
                          for name, linspace in zip(tors_names, tors_linspaces)]
            tors_sym_nums = list(automol.zmatrix.torsional_symmetry_numbers(
                zma, tors_names, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key))
            # print('tors_names:', tors_names)
            if tors_names:
                for tors_name, tors_grid, sym_num in zip(tors_names, tors_grids, tors_sym_nums):
                    locs_lst = []
                    enes = []
                    for grid_val in tors_grid:
                        locs_lst.append([[tors_name], [grid_val]])
                    for locs in locs_lst:
                        if scn_save_fs.leaf.exists(locs):
                            enes.append(scn_save_fs.leaf.file.energy.read(locs))
                        else:
                            enes.append(10.)
                            print('ERROR: missing grid value for torsional potential of {}'.
                                  format(spc_info[0]))
                    enes = numpy.subtract(enes, min_ene)
                    pot = list(enes*phycon.EH2KCAL)

                    # Build a potential list from only successful calculations
                    # print('pot test in zpe:', pot)
                    pot = _hrpot_spline_fitter(pot)

                    axis = coo_dct[tors_name][1:3]
                    atm_key = axis[1]
                    if ts_bnd:
                        for atm in axis:
                            if atm in ts_bnd:
                                atm_key = atm
                                break
                    group = list(
                        automol.graph.branch_atom_keys(gra, atm_key, axis, saddle=saddle, ts_bnd=ts_bnd) -
                        set(axis))
                    if not group:
                        for atm in axis:
                            if atm != atm_key:
                                atm_key = atm
                        group = list(
                            automol.graph.branch_atom_keys(gra, atm_key, axis, saddle=saddle, ts_bnd=ts_bnd) -
                            set(axis))
                    if saddle:
                        n_atm = automol.zmatrix.count(zma)
                        if 'addition' in spc_dct_i['class'] or 'abstraction' in spc_dct_i['class']:
                            group2 = []
                            ts_bnd1 = min(ts_bnd)
                            ts_bnd2 = max(ts_bnd)
                            for idx in range(ts_bnd2, n_atm):
                                group2.append(idx)
                            if ts_bnd1 in group:
                                for atm in group2:
                                    if atm not in group:
                                        group.append(atm)
                        # check to see if symmetry of XH3 rotor was missed
                        if sym_num == 1:
                            group2 = []
                            for idx in range(n_atm):
                                if idx not in group and idx not in axis:
                                    group2.append(idx)
                            all_H = True
                            symbols = automol.zmatrix.symbols(zma)
                            #print('symbols test:', symbols)
                            #print('second group2:', group2)
                            #print('len pot:', len(pot))
                            H_count = 0
                            for idx in group2:
                                if symbols[idx] != 'H' and symbols[idx] != 'X':
                                    all_H = False
                                    break
                                else:
                                    if symbols[idx] == 'H':
                                        H_count += 1
                            if all_H and H_count == 3:
                                sym_num = 3
                                lpot = int(len(pot)/3)
                                potp = []
                                potp[0:lpot] = pot[0:lpot]
                                pot = potp
                                #potp = []
                                #for idx in range(lpot):
                                #    potp.append(pot[idx])
                                #print('potp test:', potp)
                                #print('potp2 test:', potp)
                                #pot = potp
                            #print('all_h test=:', all_H, H_count)
                            #print('pot test:', pot)
                            #print('len pot new:', len(pot))
                    group = list(numpy.add(group, 1))
                    axis = list(numpy.add(axis, 1))
                    #print('axis test:', axis)
                    #print('atm_key:', atm_key)
                    #print('group:', group)
                    #for idx, atm in enumerate(axis):
                    #    if atm == atm_key+1:
                    #        if idx != 1:
                    #            axis.reverse()
                    #            print('axis reversed:', axis)
                    if (atm_key+1) != axis[1]:
                        axis.reverse()
                        #print('axis reversed:', axis)
                    #if atm_key != axis(0):
                        #axis.reverse()

                    #check for dummy transformations
                    atom_symbols = automol.zmatrix.symbols(zma)
                    dummy_idx = []
                    for atm_idx, atm in enumerate(atom_symbols):
                        if atm == 'X':
                            dummy_idx.append(atm_idx)
                    remdummy = numpy.zeros(len(zma[0]))
                    for dummy in dummy_idx:
                        for idx, _ in enumerate(remdummy):
                            if dummy < idx:
                                remdummy[idx] += 1

                    hind_rot_str += mess_io.writer.rotor_hindered(
                        group, axis, sym_num, pot, remdummy=remdummy)

                    #print('projrot 1 test:')
                    proj_rotors_str += projrot_io.writer.rotors(
                        axis, group, remdummy=remdummy)
                    sym_factor /= sym_num

                # Write the string for the ProjRot input
                coord_proj = 'cartesian'
                grad = ''
                #print('projrot zpe test:')
                projrot_inp_str = projrot_io.writer.rpht_input(
                    tors_geo, grad, hess, rotors_str=proj_rotors_str,
                    coord_proj=coord_proj)

                bld_locs = ['PROJROT', 0]
                bld_save_fs = autofile.fs.build(tors_save_path)
                bld_save_fs.leaf.create(bld_locs)
                path = bld_save_fs.leaf.path(bld_locs)
                print('Build Path for Partition Functions')
                print(path)
                proj_file_path = os.path.join(path, 'RPHt_input_data.dat')
                with open(proj_file_path, 'w') as proj_file:
                    proj_file.write(projrot_inp_str)

                moldr.util.run_script(projrot_script_str, path)

                zpe_har_no_tors = har_zpe
                if pot:
                    rthrproj_freqs, _ = projrot_io.reader.rpht_output(
                        path+'/hrproj_freq.dat')
                    freqs = rthrproj_freqs
                    zpe_har_no_tors = sum(freqs)*phycon.WAVEN2KCAL/2.

                # now try again with the other projrot parameters
                projrot_script_str2 = ("#!/usr/bin/env bash\n"
                "RPHt.exe >& /dev/null")
                moldr.util.run_script(projrot_script_str2, path)
                zpe_har_no_tors_2 = har_zpe
                freqs_2 = []
                if pot:
                    rthrproj_freqs_2, _ = projrot_io.reader.rpht_output(
                        path+'/hrproj_freq.dat')
                    freqs_2 = rthrproj_freqs_2
                    zpe_har_no_tors_2 = sum(freqs_2)*phycon.WAVEN2KCAL/2.

                dummy_freqs = [1000.]
                dummy_zpe = 0.0
                core = mess_io.writer.core_rigidrotor(tors_geo, sym_factor)
                spc_str = mess_io.writer.molecule(
                    core, dummy_freqs, elec_levels,
                    hind_rot=hind_rot_str,
                    )

                # create a messpf input file
                temp_step = 100.
                ntemps = 5
                zpe_str = '{0:<8.2f}\n'.format(dummy_zpe)
                zpe_str = ' ZeroEnergy[kcal/mol] ' + zpe_str
                zpe_str += 'End\n'
                global_pf_str = mess_io.writer.global_pf(
                    [], temp_step, ntemps, rel_temp_inc=0.001,
                    atom_dist_min=0.6)
                spc_head_str = 'Species ' + ' Tmp'
                pf_inp_str = '\n'.join(
                    [global_pf_str, spc_head_str,
                     spc_str, zpe_str])

                bld_locs = ['PF', 0]
                bld_save_fs = autofile.fs.build(tors_save_path)
                bld_save_fs.leaf.create(bld_locs)
                pf_path = bld_save_fs.leaf.path(bld_locs)

                # run messpf
                with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
                    pf_file.write(pf_inp_str)
                moldr.util.run_script(pf_script_str, pf_path)

                with open(os.path.join(pf_path, 'pf.log'), 'r') as mess_file:
                    output_string = mess_file.read()

                # Read the freqs and zpes
                tors_freqs = mess_io.reader.tors.freqs(output_string)
                tors_zpes = mess_io.reader.tors.zpves(output_string)
                tors_zpe_cor = 0.0
                tors_zpe = 0.0
                for (tors_freq, tors_1dhr_zpe) in zip(tors_freqs, tors_zpes):
                    tors_zpe_cor += tors_1dhr_zpe - tors_freq*phycon.WAVEN2KCAL/2
                    tors_zpe += tors_1dhr_zpe

                har_tors_zpe = har_zpe - zpe_har_no_tors
                har_tors_zpe_2 = har_zpe - zpe_har_no_tors_2
                del_tors_zpe = har_tors_zpe - tors_zpe
                del_tors_zpe_2 = har_tors_zpe_2 - tors_zpe
                #print('tors_zpe test:', del_tors_zpe, del_tors_zpe_2)
                if del_tors_zpe <= del_tors_zpe_2:
                    zpe = zpe_har_no_tors + tors_zpe
                else:
                    zpe = zpe_har_no_tors_2 + tors_zpe
                if abs(del_tors_zpe) > 0.2 and abs(del_tors_zpe_2) > 0.2:
                    print('Warning: There is a difference of {0:.2f} and {1:.2f} kcal/mol '.format(
                        del_tors_zpe, del_tors_zpe_2),
                        'between the harmonic and hindered torsional zero-point energies')
                # read torsional harmonic zpe and actual zpe
                print('zpe test in get_zpe:',zpe_har_no_tors, zpe_har_no_tors_2, 
                       tors_zpe, har_zpe, zpe)

        # used to take full harmonic zpe and add torsional hr vs harmonic diff
        #zpe = har_zpe + tors_zpe_cor
        # now take harmonic zpe for non-torsional and add tors_zpe
            #print('zpe test:', har_zpe, zpe_har_no_tors, tors_zpe, zpe)
        ret = zpe

    elif vib_model == 'HARM' and tors_model == 'MDHR':
        print('HARM and MDHR combination is not yet implemented')

    elif vib_model == 'HARM' and tors_model == 'TAU':
        print('HARM and TAU combination is not yet implemented')

    elif vib_model == 'VPT2' and tors_model == 'RIGID':
        if anh_min_cnf_locs is not None:
            anh_geo = anh_cnf_save_fs.leaf.file.geometry.read(anh_min_cnf_locs)
            min_ene = anh_cnf_save_fs.leaf.file.energy.read(anh_min_cnf_locs)
            if automol.geom.is_atom(anh_geo):
                # print('This is an atom')
                mass = ptab.to_mass(anh_geo[0][0])
                spc_str = mess_io.writer.atom(
                    mass, elec_levels)
            else:
                hess = anh_cnf_save_fs.leaf.file.hessian.read(anh_min_cnf_locs)
                freqs = elstruct.util.harmonic_frequencies(anh_geo, hess, project=True)
                if automol.geom.is_linear(anh_geo):
                    proj_freqs = freqs[5:]
                else:
                    proj_freqs = freqs[6:]

                zpe = sum(proj_freqs)*phycon.WAVEN2KCAL/2.
                hind_rot_str = ""

                core = mess_io.writer.core_rigidrotor(anh_geo, sym_factor)
                spc_str = mess_io.writer.molecule(
                    core, proj_freqs, elec_levels,
                    hind_rot=hind_rot_str,
                    )
        else:
            spc_str = ''
        print('VPT2 and RIGID combination is not yet properly implemented')

    elif vib_model == 'VPT2' and tors_model == '1DHR':
        print('VPT2 and 1DHR combination is not yet implemented')

    elif vib_model == 'VPT2' and tors_model == 'TAU':
        print('VPT2 and TAU combination is not yet implemented')

    return ret, is_atom


def tau_pf_write(
        name, save_prefix,
        run_grad=False, run_hess=False):
    """ Print out data fle for partition function evaluation
    """
    cnf_save_fs = autofile.fs.conformer(save_prefix)
    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        ene_ref = cnf_save_fs.leaf.file.energy.read(min_cnf_locs)

    tau_save_fs = autofile.fs.tau(save_prefix)
    evr = name+'\n'
    # cycle through saved tau geometries
    idx = 0
    for locs in tau_save_fs.leaf.existing():
        geo = tau_save_fs.leaf.file.geometry.read(locs)
        ene = tau_save_fs.leaf.file.energy.read(locs)
        ene = (ene - ene_ref) * phycon.EH2KCAL
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
        # print('integral convergence for T = ', temp)
        for locs in tau_save_fs.leaf.existing():
            idx += 1
            ene = tau_save_fs.leaf.file.energy.read(locs)
            ene = (ene - ene_ref) * phycon.EH2KCAL
            tmp = numpy.exp(-ene*349.7/(0.695*temp))
            sumq = sumq + tmp
            sum2 = sum2 + tmp**2
            sigma = numpy.sqrt(
                (abs(sum2/float(idx)-(sumq/float(idx))**2))/float(idx))
            print(sumq/float(idx), sigma, 100.*sigma*float(idx)/sumq, idx)


def _hrpot_spline_fitter(pot, thresh=-0.05):
    """ Get a physical hindered rotor potential via a series of spline fits
    """

    # Build a potential list from only successful calculations
    lpot = len(pot)+1
    idx_success = []
    pot_success = []
    pot.append(0.)
    for idx in range(lpot):
        if pot[idx] < 600.:
            idx_success.append(idx)
            pot_success.append(pot[idx])
    idx_success.append(lpot)
    pot_success.append(pot[0])
    pot_spl = interp1d(
        numpy.array(idx_success), numpy.array(pot_success), kind='cubic')
    for idx in range(lpot):
        pot[idx] = pot_spl(idx)

    # Do second spline fit of only positive values if any negative values found
    if any(val < thresh for val in pot):
        print('Found pot vals below {0} kcal. Refit w/ positives'.format(thresh))
        print('Potential before spline:', pot)
        x_pos = numpy.array([i for i in range(lpot)
                             if pot[i] >= thresh])
        y_pos = numpy.array([pot[i] for i in range(lpot)
                             if pot[i] >= thresh])
        pos_pot_spl = interp1d(x_pos, y_pos, kind='cubic')
        pot_pos_fit = []
        for idx in range(lpot):
            pot_pos_fit.append(pos_pot_spl(idx))

        print('Potential after spline:', pot_pos_fit)
        # Perform second check to see if negative potentials have been fixed
        if any(val < thresh for val in pot_pos_fit):
            print('Found values below {0} kcal again. Trying linear interp of positive vals'
                  .format(thresh))
            neg_idxs = [i for i in range(lpot) if pot_pos_fit[i] < thresh]
            clean_pot = []
            for i in range(lpot):
                if i in neg_idxs:
                    # Find the indices for positive vals around negative value
                    idx_0 = i - 1
                    while idx_0 in neg_idxs:
                        idx_0 = idx_0 - 1
                    for j in range(i, lpot):
                        if pot_pos_fit[j] >= thresh:
                            idx_1 = j
                            break
                    # Get a new value for this point on the potential by
                    # doing a linear interp of positives
                    interp_val = (
                        pot_pos_fit[idx_0] * (1.0 - ((i - idx_0) / (idx_1 - idx_0))) +
                        pot_pos_fit[idx_1] * ((i - idx_0) / (idx_1 - idx_0))
                    )
                    clean_pot.append(interp_val)
                else:
                    clean_pot.append(pot[i])
            final_potential = clean_pot.copy()

        else:
            final_potential = pot_pos_fit.copy()

    else:
        final_potential = pot.copy()

    print('Final potential in spline fittere:', final_potential)
    final_potential = final_potential[:-1]

    return final_potential
