""" drivers
"""
import os
import numpy
from qcelemental import constants as qcc
from qcelemental import periodictable as ptab
import projrot_io
import automol
import elstruct
import autofile
import moldr
import mess_io

WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

def species_block(
        spc, spc_dct_i, spc_info, spc_model, pf_levels, projrot_script_str,
        elec_levels=[[0., 1]], sym_factor=1.,
        save_prefix='spc_save_path'):
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
    if 'ts_' in spc:
        har_save_fs = autofile.fs.ts(har_save_path)
        har_save_fs.trunk.create()
        har_save_path = har_save_fs.trunk.path()
        saddle = True
        if 'migration' in spc_dct_i[spc]['class'] or 'elimination' in spc_dct_i[spc]['class']:
            dist_names.append(spc_dct_i[spc]['dist_info'][0])
            dist_names.append(spc_dct_i[spc]['dist_info'][3])
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

    if tors_level:
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
                form_coords = list(automol.bond_idxs(zma, dist_names[0]))
                form_coords.extend(list(automol.bond_idxs(zma, dist_names[1])))
            sym_factor = moldr.conformer.symmetry_factor(sym_geo, sym_ene, sym_cnf_save_fs, saddle, form_coords)
            # xyzs = automol.geom.coordinates(sym_geo)
            print('sym_factor from moldr sampling:', sym_factor)
        if sym_model == '1DHR':
            # Warning: the 1DHR based symmetry number has not yet been set up
            sym_factor = 1

    imag_freq = 0.
    print('sym_factor calculated:', sym_factor)

    if vib_model == 'HARM' and tors_model == 'RIGID':
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

    if vib_model == 'HARM' and tors_model == '1DHR':
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
                hind_rot_str = ""
                proj_rotors_str = ""

                if tors_min_cnf_locs is not None:
                    if har_cnf_save_fs.trunk.file.info.exists():
                        inf_obj_s = har_cnf_save_fs.trunk.file.info.read()
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
                        ts_bnd = automol.zmatrix.bond_idxs(zma, dist_name)
                        ts_bnd = frozenset(ts_bnd)
                    pot = []
                    if 'hind_inc' in spc_dct_i:
                        scan_increment = spc_dct_i['hind_inc']
                    else:
                        scan_increment = 30. * qcc.conversion_factor('degree', 'radian')
                    val_dct = automol.zmatrix.values(zma)
                    tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                        zma, tors_names, scan_increment)
                    tors_grids = [
                        numpy.linspace(*linspace) + val_dct[name]
                        for name, linspace in zip(tors_names, tors_linspaces)]
                    tors_sym_nums = list(automol.zmatrix.torsional_symmetry_numbers(zma, tors_names))
                    print('tors sym nums:', tors_sym_nums)
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
                                enes.append(20.)
                                print('ERROR: missing grid value for torsional potential of {}'.format(spc_info[0]))
                        enes = numpy.subtract(enes, min_ene)
                        pot = list(enes*EH2KCAL)
                        print('pot test in species_block:', pot)
                        axis = coo_dct[tors_name][1:3]
                        group = list(
                            automol.graph.branch_atom_keys(gra, axis[1], axis, saddle=saddle, ts_bnd=ts_bnd) -
                            set(axis))
                        group = list(numpy.add(group, 1))
                        axis = list(numpy.add(axis, 1))

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
                            group, axis, sym_num, pot, remdummy)
                        proj_rotors_str += projrot_io.writer.rotors(
                            axis, group, remdummy=remdummy)
                        sym_factor /= sym_num
                        print('sym_factor renorm:', sym_factor)
                        idx += 1

                    # Write the string for the ProjRot input
                    coord_proj = 'cartesian'
                    grad = ''
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

                    freqs = []
                    if len(pot) > 0:
                        rthrproj_freqs, imag_freq = projrot_io.reader.rpht_output(
                            path+'/hrproj_freq.dat')
                        freqs = rthrproj_freqs
                    if not freqs:
                        rtproj_freqs, imag_freq = projrot_io.reader.rpht_output(
                            path+'/RTproj_freq.dat')
                        freqs = rtproj_freqs
                    if 'ts_' in spc:
                        if imag_freq:
                           imag_freq = imag_freq[0]
                        else:
                            imag_freq = freqs[-1]
                            freqs = freqs[:-1]
                    # shutil.rmtree(path)
                    core = mess_io.writer.core_rigidrotor(tors_geo, sym_factor)
                    spc_str = mess_io.writer.molecule(
                        core, freqs, elec_levels,
                        hind_rot=hind_rot_str
                        )
        else:
            print('ERROR: Reference geometry is missing for harmonic frequencies for species {}'.format(spc_info[0]))
            spc_str = ''
            imag_freq = 0.
            return '', 0.

    if vib_model == 'HARM' and tors_model == 'MDHR':
        print('HARM and MDHR combination is not yet implemented')

    if vib_model == 'HARM' and tors_model == 'TAU':
        print('HARM and TAU combination is not yet implemented')
        moldr.driver.tau_pf_write(
            name=name,
            save_prefix=thy_save_path,
            run_grad=run_grad_pf,
            run_hess=run_hess_pf,
        )
    if vib_model == 'VPT2' and tors_model == 'RIGID':
        if anh_min_cnf_locs is not None:
            anh_geo = anh_cnf_save_fs.leaf.file.geometry.read(anh_min_cnf_locs)
            min_ene = anh_cnf_save_fs.leaf.file.energy.read(anh_min_cnf_locs)
            if automol.geom.is_atom(anh_geo):
                print('This is an atom')
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
                zpe = sum(freqs)*WAVEN2KCAL/2.
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

    if vib_model == 'VPT2' and tors_model == '1DHR':
        print('VPT2 and 1DHR combination is not yet implemented')

    if vib_model == 'VPT2' and tors_model == 'TAU':
        print('VPT2 and TAU combination is not yet implemented')
        moldr.driver.tau_pf_write(
            name=name,
            save_prefix=thy_save_path,
            run_grad=run_grad_pf,
            run_hess=run_hess_pf,
        )

    return spc_str, imag_freq


def vtst_with_saddle_block(
        ts_dct,  ts_label, rct_label, prd_label, spc_ene, rct_zpe, projrot_script_str,
        multi_info, elec_levels=[[0., 1]], sym_factor=1.
        ):
    """ prepare the mess input string for a variational TS that has
    a saddle point. Do it by calling the species block for each grid point
    in the scan file system
    """

    orb_restr = moldr.util.orbital_restriction(ts_info, multi_info)
    multi_level = multi_info[0:3]
    multi_level.append(orb_restr)

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

    for idx, grid_val in enumerate(grid):
        locs = [[dist_name], [grid_val]]
        if not scn_save_fs.leaf.file.geometry.exists(locs):
            continue
        else:
            geom = scn_save_fs.leaf.file.geometry.read(locs)
        if not scn_save_fs.leaf.file.energy.exists(locs):
            continue
        else:
            ene = scn_save_fs.leaf.file.energy.read(locs)
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
                rthrproj_freqs, imag_freq = projrot_io.reader.rpht_output(
                    path+'/hrproj_freq.dat')
                freqs = rthrproj_freqs
            if not freqs:
                rtproj_freqs, imag_freq = projrot_io.reader.rpht_output(
                    path+'/RTproj_freq.dat')
                freqs = rtproj_freqs
                if not imag_freq:
                    freqs = freqs[:-1]

        zpe = sum(freqs)*WAVEN2KCAL/2.
        erel = (ene - inf_sep_ene)*EH2KCAL
        erel_zpe_corr = erel + zpe - rct_zpe
        eref = erel - spc_ene

        # Iniialize the header of the string
        irc_pt_str = '!-----------------------------------------------'
        irc_pt_str += '! IRC Point {0}\n'.format(str(idx+1))

        # Write the molecule section for each irc point
        core = mess_io.writer.mol_data.core_rigidrotor(geom, sym_factor, interp_emax='')
        irc_pt_str += mess_io.writer.species.molecule(core, freqs, elec_levels,
             hind_rot='', xmat=None, rovib_coups='', rot_dists='')

        # Append the zero point energy for the molecule
        irc_pt_str += '    ZeroEnergy[kcal/mol]      {0:}'.format(eref)

        # Append string to list
        irc_pt_strs.append(irc_pt_str)

    # Write the MESS string for the variational sections
    varational_str = mess_io.writer.rxnchan.ts_variational(
        ts_label, reac_label, prod_label, irc_pt_strs)
    print('variational_str test:', variational_str)

    return variational_str


def vtst_no_saddle_block(scn_save_fs, geoms, frequencies, energies):
    """ prepare the mess input string for a variational TS where there is not a 
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
        irc_pt_str += '    ZeroEnergy[kcal/mol]      {0:}'.format(zero_energy)

        # Append string to list
        irc_pt_strings.append(irc_pt_string)

    # Write the MESS string for the variational sections
    varational_str = mess_io.writer.rxnchan.ts_variational(
        ts_label, reac_label, prod_label, irc_pt_strings)

    return variational_str


def pst_block(
        spc_dct_i, spc_dct_j, spc_model, pf_levels, projrot_script_str,
        spc_save_fs, elec_levels=[[0., 1]], sym_factor=1.
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
        elec_levs_i = spc_dct_i['elec_levs']
    else:
        elec_levs_i = [[0., spc_dct_i['mul']]]
    if 'elec_levs' in spc_dct_j:
        elec_levs_j = spc_dct_j['elec_levs']
    else:
        elec_levs_j = [[0., spc_dct_j['mul']]]

    elec_levs = []
    for _, elec_lev_i in enumerate(elec_levs_i):
        for _, elec_lev_j in enumerate(elec_levs_j):
            elec_levs.append(
                [elec_lev_i[0]+elec_lev_j[0],
                 elec_lev_i[1]*elec_lev_j[1]])

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
            sym_factor_ = 1
    sym_factor = sym_factor_i * sym_factor_j

    if vib_model:
    # if vib_model == 'HARM' and tors_model == 'RIGID':
        if har_min_cnf_locs_i is not None:
            har_geo_i = har_cnf_save_fs_i.leaf.file.geometry.read(har_min_cnf_locs_i)
            if har_min_cnf_locs_j is not None:
                har_geo_j = har_cnf_save_fs_j.leaf.file.geometry.read(har_min_cnf_locs_j)
                if not automol.geom.is_atom(har_geo_i):
                    hess_i = har_cnf_save_fs_i.leaf.file.hessian.read(har_min_cnf_locs_i)
                    freqs_i = elstruct.util.harmonic_frequencies(har_geo_i, hess_i, project=False)
                    mode_start = 6
                    if automol.geom.is_linear(har_geo_i):
                        mode_start = mode_start - 1
                    freqs = freqs_i[mode_start:]
                    if not automol.geom.is_atom(har_geo_j):
                        hess_j = har_cnf_save_fs_j.leaf.file.hessian.read(har_min_cnf_locs_j)
                        freqs_j = elstruct.util.harmonic_frequencies(har_geo_j, hess_j, project=False)
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
                    har_geo_i, har_geo_j, sym_factor, stoich, pot_prefactor=10., pot_power_exp=6)
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


    spc_str = ''
        # elec_levels_j = spc_dct_j['elec_levs']

    elec_levels = [[0., 2]]
    sym_factor = 1.

    if vib_model == 'HARM' and tors_model == 'RIGID':
        if har_min_cnf_locs_i is not None:
            har_geo_i = har_cnf_save_fs_i.leaf.file.geometry.read(har_min_cnf_locs_i)
            if har_min_cnf_locs_j is not None:
                har_geo_j = har_cnf_save_fs_j.leaf.file.geometry.read(har_min_cnf_locs_j)

                freqs = (30, 50, 70, 100, 200)
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

                har_geo = har_geo_i
                har_geo_j = automol.geom.translated(har_geo_j, [10., 10., 10.])
                har_geo += har_geo_j

                hind_rot_str = ""

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

    har_cnf_save_fs = autofile.fs.conformer(har_save_path)
    har_min_cnf_locs = moldr.util.min_energy_conformer_locators(har_cnf_save_fs)

    if tors_level:
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
        full_freqs = elstruct.util.harmonic_frequencies(har_geo, hess, project=False)

        mode_start = 6
        if 'ts_' in spc:
            mode_start = mode_start + 1
        if automol.geom.is_linear(har_geo):
            mode_start = mode_start - 1
        proj_freqs = full_freqs[mode_start:]

        har_zpe = sum(proj_freqs)*WAVEN2KCAL/2.

    if vib_model == 'HARM' and tors_model == 'RIGID':
        ret = har_zpe

    if vib_model == 'HARM' and tors_model == '1DHR':
        # make pf string for 1d rotor
        # run messpf
        # read 1d harmonic and torsional ZPEs
        # modify har_zpe

        hind_rot_str = ""
        if tors_min_cnf_locs is not None:
            if har_cnf_save_fs.trunk.file.info.exists():
                inf_obj_s = har_cnf_save_fs.trunk.file.info.read()
                tors_ranges = inf_obj_s.tors_ranges
                tors_ranges = autofile.info.dict_(tors_ranges)
                print(tors_ranges)
                tors_names = list(tors_ranges.keys())
            else:
                print('No inf obj to identify torsional angles')
                tors_names = []

        min_ene = tors_cnf_save_fs.leaf.file.energy.read(tors_min_cnf_locs)
        tors_geo = tors_cnf_save_fs.leaf.file.geometry.read(tors_min_cnf_locs)
        zma = tors_cnf_save_fs.leaf.file.zmatrix.read(tors_min_cnf_locs)
        gra = automol.zmatrix.graph(zma, remove_stereo=True)
        tors_zpe_cor = 0.0
        if tors_names:
            coo_dct = automol.zmatrix.coordinates(zma, multi=False)
            # prepare axis, group, info
            scn_save_fs = autofile.fs.scan(tors_cnf_save_path)
            ts_bnd = None
            if saddle:
                dist_name = spc_dct_i['dist_info'][0]
                ts_bnd =  automol.zmatrix.bond_idxs(zma, dist_name)
                ts_bnd = frozenset(ts_bnd)
            pot = []
            if 'hind_inc' in spc_dct_i:
                scan_increment = spc_dct_i['hind_inc']
            else:
                scan_increment = 30. * qcc.conversion_factor('degree', 'radian')
            val_dct = automol.zmatrix.values(zma)
            tors_linspaces = automol.zmatrix.torsional_scan_linspaces(
                zma, tors_names, scan_increment)
            tors_grids = [numpy.linspace(*linspace) + val_dct[name]
                          for name, linspace in zip(tors_names, tors_linspaces)]
            tors_sym_nums = list(automol.zmatrix.torsional_symmetry_numbers(zma, tors_names))
            for tors_name, tors_grid, sym_num in zip(tors_names, tors_grids, tors_sym_nums):
                print('tors test:', tors_name, tors_grid, sym_num)
                locs_lst = []
                enes = []
                for grid_val in tors_grid:
                    locs_lst.append([[tors_name], [grid_val]])
                for locs in locs_lst:
                    if scn_save_fs.leaf.exists(locs):
                        enes.append(scn_save_fs.leaf.file.energy.read(locs))
                    else:
                        enes.append(1000.)
                        print('ERROR: missing grid value for torionsal potential of {}'.format(spc_info[0]))
                enes = numpy.subtract(enes, min_ene)
                pot = list(enes*EH2KCAL)
                print('pot test in zero point energy:', pot)
                axis = coo_dct[tors_name][1:3]
                group = list(
                    automol.graph.branch_atom_keys(gra, axis[1], axis, saddle=saddle, ts_bnd=ts_bnd) -
                    set(axis))
                group = list(numpy.add(group, 1))
                axis = list(numpy.add(axis, 1))

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
                    group, axis, sym_num, pot, remdummy)

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
            for (tors_freq, tors_1dhr_zpe) in zip(tors_freqs, tors_zpes):
                tors_zpe_cor += tors_1dhr_zpe - tors_freq*WAVEN2KCAL/2

            # read torsional harmonic zpe and actual zpe

        zpe = har_zpe + tors_zpe_cor
        ret = zpe

    if vib_model == 'HARM' and tors_model == 'MDHR':
        print('HARM and MDHR combination is not yet implemented')

    if vib_model == 'HARM' and tors_model == 'TAU':
        print('HARM and TAU combination is not yet implemented')

    if vib_model == 'VPT2' and tors_model == 'RIGID':
        if anh_min_cnf_locs is not None:
            anh_geo = anh_cnf_save_fs.leaf.file.geometry.read(anh_min_cnf_locs)
            min_ene = anh_cnf_save_fs.leaf.file.energy.read(anh_min_cnf_locs)
            if automol.geom.is_atom(anh_geo):
                print('This is an atom')
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

                zpe = sum(proj_freqs)*WAVEN2KCAL/2.
                hind_rot_str = ""

                core = mess_io.writer.core_rigidrotor(anh_geo, sym_factor)
                spc_str = mess_io.writer.molecule(
                    core, proj_freqs, elec_levels,
                    hind_rot=hind_rot_str,
                    )
        else:
            spc_str = ''
        print('VPT2 and RIGID combination is not yet properly implemented')

    if vib_model == 'VPT2' and tors_model == '1DHR':
        print('VPT2 and 1DHR combination is not yet implemented')

    if vib_model == 'VPT2' and tors_model == 'TAU':
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
