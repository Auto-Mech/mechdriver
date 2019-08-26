""" drivers
"""
import os
import functools
import numpy
from qcelemental import constants as qcc
from qcelemental import periodictable as ptab
import projrot_io
import automol
import elstruct
import autofile
import moldr
from moldr import runner
import mess_io

WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

def species_block(
        spc_info,
        tors_model, vib_model, 
        har_level, tors_level, vpt2_level,
        script_str,
#        orb_restr, script_str,
        elec_levels=[[0., 1]], sym_factor=1.,
        save_prefix='spc_save_path'):
    """ prepare the species input for messpf
    """

    # prepare the three sets of file systems
    orb_restr = moldr.util.orbital_restriction(
        spc_info, har_level)
    har_levelp = har_level[1:3]
    har_levelp.append(orb_restr)

    thy_save_fs = autofile.fs.theory(save_prefix)

    har_save_path = thy_save_fs.leaf.path(har_levelp)
    print('har_save_path:',har_save_path)
    har_min_cnf_locs = moldr.util.min_energy_conformer_locators(har_save_path)
    har_cnf_save_fs = autofile.fs.conformer(har_save_path)
    har_cnf_save_path = har_cnf_save_fs.leaf.path(har_min_cnf_locs)

    orb_restr = moldr.util.orbital_restriction(
        spc_info, tors_level)
    tors_levelp = tors_level[1:3]
    tors_levelp.append(orb_restr)

    tors_save_path = thy_save_fs.leaf.path(tors_levelp)
    tors_min_cnf_locs = moldr.util.min_energy_conformer_locators(tors_save_path)
    tors_cnf_save_fs = autofile.fs.conformer(tors_save_path)
    tors_cnf_save_path = tors_cnf_save_fs.leaf.path(tors_min_cnf_locs)

    orb_restr = moldr.util.orbital_restriction(
        spc_info, vpt2_level)
    vpt2_levelp = vpt2_level[1:3]
    vpt2_levelp.append(orb_restr)

    anh_save_path = thy_save_fs.leaf.path(vpt2_levelp)
    anh_min_cnf_locs = moldr.util.min_energy_conformer_locators(anh_save_path)
    anh_cnf_save_fs = autofile.fs.conformer(anh_save_path)
    anh_cnf_save_path = anh_cnf_save_fs.leaf.path(anh_min_cnf_locs)

    # atom case - do as first step in each of other cases
    # pure harmonic case
    spc_str = ''
    if vib_model == 'HARM' and tors_model == 'RIGID':
        if har_min_cnf_locs is not None:
            har_geo = har_cnf_save_fs.leaf.file.geometry.read(har_min_cnf_locs)
            min_ene = har_cnf_save_fs.leaf.file.energy.read(har_min_cnf_locs)
            if automol.geom.is_atom(har_geo):
                print('This is an atom')
                mass = ptab.to_mass(har_geo[0][0])
                spc_str = mess_io.writer.write_atom(
                    mass, elec_levels)
            else:
                grad = har_cnf_save_fs.leaf.file.gradient.read(har_min_cnf_locs)
                hess = har_cnf_save_fs.leaf.file.hessian.read(har_min_cnf_locs)
                freqs = elstruct.util.harmonic_frequencies(har_geo, hess, project=False)
#                freqs = elstruct.util.harmonic_frequencies(har_geo, hess, project=True)
                if automol.geom.is_linear(har_geo):
                    proj_freqs = freqs[5:]
                else:
                    proj_freqs = freqs[6:]

                print('projected freqs including low frequencies')
                print(freqs)
                print('projected freqs')
                print(proj_freqs)
                zpe = sum(proj_freqs)*WAVEN2KCAL/2.
                hind_rot_str = ""

                core = mess_io.writer.write_core_rigidrotor(har_geo, sym_factor)
                spc_str = mess_io.writer.write_molecule(
                    core, proj_freqs, elec_levels,
                    hind_rot=hind_rot_str,
                    )
        else:
            spc_str = ''

    if vib_model == 'HARM' and tors_model == '1DHR':
        if har_min_cnf_locs is not None:
            har_geo = har_cnf_save_fs.leaf.file.geometry.read(har_min_cnf_locs)
            min_ene = har_cnf_save_fs.leaf.file.energy.read(har_min_cnf_locs)
            if automol.geom.is_atom(har_geo):
                print('This is an atom')
                mass = ptab.to_mass(har_geo[0][0])
                spc_str = mess_io.writer.write_atom(
                    mass, elec_levels)
            else:
                grad = har_cnf_save_fs.leaf.file.gradient.read(har_min_cnf_locs)
                hess = har_cnf_save_fs.leaf.file.hessian.read(har_min_cnf_locs)
                freqs = elstruct.util.harmonic_frequencies(har_geo, hess, project=False)
                # freqs = elstruct.util.harmonic_frequencies(har_geo, hess, project=True)
                print(automol.geom.is_linear(har_geo))
                hind_rot_str = ""
                proj_rotors_str = ""

                if tors_min_cnf_locs is not None:
                    tors_geo = tors_cnf_save_fs.leaf.file.geometry.read(tors_min_cnf_locs)
                    if automol.geom.is_linear(har_geo):
                        proj_freqs = freqs[5:]
                        zpe = sum(proj_freqs)*WAVEN2KCAL/2.
                    else:
                        zma = automol.geom.zmatrix(tors_geo)
                        gra = automol.zmatrix.graph(zma, remove_stereo=True)
                        tors_names = automol.geom.zmatrix_torsion_coordinate_names(tors_geo)
                        coo_dct = automol.zmatrix.coordinates(zma, multi=False)

                        # prepare axis, group, and projection info
                        scn_save_fs = autofile.fs.scan(tors_cnf_save_path)
                        print('tors_names')
                        print(tors_names)
                        pot = []
                        for tors_name in tors_names:
                            enes = [scn_save_fs.leaf.file.energy.read(locs) for locs
                                    in scn_save_fs.leaf.existing([[tors_name]])]
                            enes = numpy.subtract(enes, min_ene)
                            pot = list(enes*EH2KCAL)
                            axis = coo_dct[tors_name][1:3]
                            group = list(
                                automol.graph.branch_atom_keys(gra, axis[0], axis) -
                                set(axis))
                            group = list(numpy.add(group, 1))
                            axis = list(numpy.add(axis, 1))
                            sym = 1
                            hind_rot_str += mess_io.writer.write_rotor_hindered(
                                group, axis, sym, pot)
                            proj_rotors_str += projrot_io._write.write_rotors_str(
                                axis, group)
                        print('hind_rot_str test')
                        print(hind_rot_str)
                        print(pot)

                        # Write the string for the ProjRot input
                        COORD_PROJ = 'cartesian'
                        projrot_inp_str = projrot_io._write.write_rpht_input(
                            tors_geo, grad, hess, rotors_str=proj_rotors_str,
                            coord_proj=COORD_PROJ)

                        bld_locs = ['PROJROT', 0]
                        bld_save_fs = autofile.fs.build(tors_save_path)
                        bld_save_fs.leaf.create(bld_locs)
                        path = bld_save_fs.leaf.path(bld_locs)
                        print('Build Path for Partition Functions')
                        print(path)
                        proj_file_path = os.path.join(path, 'RPHt_input_data.dat')
                        with open(proj_file_path, 'w') as proj_file:
                            proj_file.write(projrot_inp_str)

                        moldr.util.run_script(script_str, path)

                        rtproj_freqs, _ = projrot_io._read.read_rpht_output(
                            path+'/RTproj_freq.dat')
                        rthrproj_freqs, _ = projrot_io._read.read_rpht_output(
                            path+'/hrproj_freq.dat')
                        # the second variable above is the imaginary frequency list
                        print('Projection test')
                        print(rtproj_freqs)
                        print(rthrproj_freqs)
                        if pot is None:
                            proj_freqs = rtproj_freqs
                            zpe = sum(rtproj_freqs)*WAVEN2KCAL/2.
                        else:
                            proj_freqs = rthrproj_freqs
                            zpe = sum(rthrproj_freqs)*WAVEN2KCAL/2.

                    core = mess_io.writer.write_core_rigidrotor(tors_geo, sym_factor)
                    spc_str = mess_io.writer.write_molecule(
                        core, proj_freqs, elec_levels,
                        hind_rot=hind_rot_str,
                        )

        else:
            spc_str = ''

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
                spc_str = mess_io.writer.write_atom(
                    mass, elec_levels)
            else:
                hess = anh_cnf_save_fs.leaf.file.hessian.read(anh_min_cnf_locs)
                freqs = elstruct.util.harmonic_frequencies(anh_geo, hess, project=True)
                if automol.geom.is_linear(anh_geo):
                    proj_freqs = freqs[5:]
                else:
                    proj_freqs = freqs[6:]

                print('projected freqs including low frequencies')
                print(freqs)
                print('projected freqs')
                print(proj_freqs)
                zpe = sum(proj_freqs)*WAVEN2KCAL/2.
                hind_rot_str = ""

                core = mess_io.writer.write_core_rigidrotor(anh_geo, sym_factor)
                spc_str = mess_io.writer.write_molecule(
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
        moldr.driver.tau_pf_write(
            name=name,
            save_prefix=thy_save_path,
            run_grad=run_grad_pf,
            run_hess=run_hess_pf,
        )

    return spc_str


def get_high_level_energy(
        spc_info, theory_low_level, theory_high_level, save_prefix):
    """ get high level energy at low level optimized geometry
    """

    spc_save_fs = autofile.fs.species(save_prefix)
    spc_save_fs.leaf.create(spc_info)
    spc_save_path = spc_save_fs.leaf.path(spc_info)

    orb_restr = moldr.util.orbital_restriction(
        spc_info, theory_low_level)
    thy_low_level = theory_low_level[1:3]
    thy_low_level.append(orb_restr)

    ll_save_fs = autofile.fs.theory(spc_save_path)
    ll_save_fs.leaf.create(thy_low_level)
    ll_save_path = ll_save_fs.leaf.path(thy_low_level)

    min_cnf_locs = moldr.util.min_energy_conformer_locators(
        ll_save_path)
    cnf_save_fs = autofile.fs.conformer(ll_save_path)
    cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
    min_cnf_geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)

    orb_restr = moldr.util.orbital_restriction(
        spc_info, theory_high_level)
    thy_high_level = theory_high_level[1:3]
    thy_high_level.append(orb_restr)

    sp_save_fs = autofile.fs.single_point(cnf_save_path)
    sp_save_fs.leaf.create(thy_high_level)

    min_ene = sp_save_fs.leaf.file.energy.read(thy_high_level)
    print('high level energy test')
    print(min_ene)

    return min_ene


def get_zero_point_energy(
        spc_info, tors_model, vib_model,
        har_level, tors_level, vpt2_level,
        script_str,
        elec_levels=[[0., 1]], sym_factor=1.,
        save_prefix='spc_save_path'):
    """ compute the ZPE including torsional and anharmonic corrections
    """

    # prepare the three sets of file systems
    print('save_prefix in get zpe')
    print(save_prefix)
    thy_save_fs = autofile.fs.theory(save_prefix)

    orb_restr = moldr.util.orbital_restriction(
        spc_info, har_level)
    har_levelp = har_level[1:3]
    har_levelp.append(orb_restr)

    har_save_path = thy_save_fs.leaf.path(har_levelp)
    har_min_cnf_locs = moldr.util.min_energy_conformer_locators(har_save_path)
    har_cnf_save_fs = autofile.fs.conformer(har_save_path)

    orb_restr = moldr.util.orbital_restriction(
        spc_info, tors_level)
    tors_levelp = tors_level[1:3]
    tors_levelp.append(orb_restr)

    tors_save_path = thy_save_fs.leaf.path(tors_levelp)
    tors_min_cnf_locs = moldr.util.min_energy_conformer_locators(tors_save_path)
    tors_cnf_save_fs = autofile.fs.conformer(tors_save_path)
    tors_cnf_save_path = tors_cnf_save_fs.leaf.path(tors_min_cnf_locs)

    orb_restr = moldr.util.orbital_restriction(
        spc_info, vpt2_level)
    vpt2_levelp = vpt2_level[1:3]
    vpt2_levelp.append(orb_restr)

    anh_save_path = thy_save_fs.leaf.path(vpt2_levelp)
    anh_min_cnf_locs = moldr.util.min_energy_conformer_locators(anh_save_path)
    anh_cnf_save_fs = autofile.fs.conformer(anh_save_path)
    anh_cnf_save_path = anh_cnf_save_fs.leaf.path(anh_min_cnf_locs)

    har_zpe = 0.0
    is_atom = False
    # get reference harmonic
    har_geo = har_cnf_save_fs.leaf.file.geometry.read(har_min_cnf_locs)
    if automol.geom.is_atom(har_geo):
        har_zpe = 0.0
        is_atom = True

    else:
        hess = har_cnf_save_fs.leaf.file.hessian.read(har_min_cnf_locs)
        full_freqs = elstruct.util.harmonic_frequencies(har_geo, hess, project=False)
        freqs = elstruct.util.harmonic_frequencies(har_geo, hess, project=True)
        if automol.geom.is_linear(har_geo):
            proj_freqs = full_freqs[5:]
            # proj_freqs = freqs[5:]
        else:
            proj_freqs = full_freqs[6:]
            # proj_freqs = freqs[6:]
        har_zpe = sum(proj_freqs)*WAVEN2KCAL/2.
        print('har zpe test')
        print(har_zpe)

    if vib_model == 'HARM' and tors_model == 'RIGID':
        ret = har_zpe

    if vib_model == 'HARM' and tors_model == '1DHR':
        # make pf string for 1d rotor
        # run messpf
        # read 1d harmonic and torsional ZPEs
        # modify har_zpe

        hind_rot_str = ""

        if tors_min_cnf_locs is not None:
            min_ene = tors_cnf_save_fs.leaf.file.energy.read(tors_min_cnf_locs)
            tors_geo = tors_cnf_save_fs.leaf.file.geometry.read(tors_min_cnf_locs)
            zma = automol.geom.zmatrix(tors_geo)
            gra = automol.zmatrix.graph(zma, remove_stereo=True)
            tors_names = automol.geom.zmatrix_torsion_coordinate_names(tors_geo)
            tors_zpe_cor = 0.0
            if tors_names is not None:
                coo_dct = automol.zmatrix.coordinates(zma, multi=False)

                # prepare axis, group, info
                scn_save_fs = autofile.fs.scan(tors_cnf_save_path)
                pot = []
                for tors_name in tors_names:
                    enes = [scn_save_fs.leaf.file.energy.read(locs) for locs
                            in scn_save_fs.leaf.existing([[tors_name]])]
                    enes = numpy.subtract(enes, min_ene)
                    pot = list(enes*EH2KCAL)
                    axis = coo_dct[tors_name][1:3]
                    group = list(
                        automol.graph.branch_atom_keys(gra, axis[0], axis) -
                        set(axis))
                    group = list(numpy.add(group, 1))
                    axis = list(numpy.add(axis, 1))
                    sym = 1
                    hind_rot_str += mess_io.writer.write_rotor_hindered(
                        group, axis, sym, pot)

                dummy_freqs = [1000.]
                dummy_zpe = 0.0
                core = mess_io.writer.write_core_rigidrotor(tors_geo, sym_factor)
                print('mess writer in get zpe')
                print(core)
                print(elec_levels)
                print(hind_rot_str)
                spc_str = mess_io.writer.write_molecule(
                    core, dummy_freqs, elec_levels,
                    hind_rot=hind_rot_str,
                    )

                # create a messpf input file
                temp_step = 100.
                ntemps = 5
                global_pf_str = mess_io.writer.write_global_pf(
                    [], temp_step, ntemps, rel_temp_inc=0.001,
                    atom_dist_min=0.6)
                spc_head_str = 'Species ' + ' Tmp'
                pf_inp_str = '\n'.join(
                    [global_pf_str, spc_head_str,
                     spc_str])

                bld_locs = ['PF', 0]
                bld_save_fs = autofile.fs.build(tors_save_path)
                bld_save_fs.leaf.create(bld_locs)
                pf_path = bld_save_fs.leaf.path(bld_locs)

                # run messpf
                with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
                    pf_file.write(pf_inp_str)
                moldr.util.run_script(script_str, pf_path)

                with open(os.path.join(pf_path, 'pf.log'), 'r') as mess_file:
                    output_string = mess_file.read()

                # Read the freqs and zpes
                tors_freqs = mess_io.reader.tors.read_freqs(output_string)
                tors_zpes = mess_io.reader.tors.read_zpes(output_string)
                tors_zpe_cor = 0.0
                print('tors zpe test')
                for (tors_freq, tors_1dhr_zpe) in zip(tors_freqs, tors_zpes):
                    tors_zpe_cor += tors_1dhr_zpe - tors_freq*WAVEN2KCAL/2
                    print(tors_1dhr_zpe, tors_freq, tors_freq*WAVEN2KCAL/2)

                # read torsional harmonic zpe and actual zpe

            zpe = har_zpe + tors_zpe_cor
            print (zpe,har_zpe,tors_zpe_cor)
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
                spc_str = mess_io.writer.write_atom(
                    mass, elec_levels)
            else:
                hess = anh_cnf_save_fs.leaf.file.hessian.read(anh_min_cnf_locs)
                freqs = elstruct.util.harmonic_frequencies(anh_geo, hess, project=True)
                if automol.geom.is_linear(anh_geo):
                    proj_freqs = freqs[5:]
                else:
                    proj_freqs = freqs[6:]

                print('projected freqs including low frequencies')
                print(freqs)
                print('projected freqs')
                print(proj_freqs)
                zpe = sum(proj_freqs)*WAVEN2KCAL/2.
                hind_rot_str = ""

                core = mess_io.writer.write_core_rigidrotor(anh_geo, sym_factor)
                spc_str = mess_io.writer.write_molecule(
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


