""" reaction list test
"""
import os
import sys
from qcelemental import constants as qcc
import thermo
import automol
import autofile
import moldr
import mess_io.writer
import scripts

ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

def species_thermo(
        spc_names,
        spc_info,
        spc_ref_names,
        elc_deg_dct,
        temp_par,
        ref_high_level,
        run_high_levels,
        spc_models,
        pf_levels,
        save_prefix,
        projrot_script_str,
        pf_script_str,
        ):
    """Evaluate the thermodynamics properties for all the species
    """

    # initializations
    temp_step = temp_par[0]
    ntemps = temp_par[1]
    chemkin_poly_strs = {}
    for hl_idx, _ in enumerate(run_high_levels):
        chemkin_poly_strs[hl_idx] = ''
    spc_ref_ich = []
    for ref_name in spc_ref_names:
        spc_ref_ich.append(spc_info[ref_name][0])

    # generate and store the energies so that the reference energies can be used as needed
    ene_hl = {}
    print('Evaluating thermo for following species')
    print(spc_names)
    for name in spc_names:
        # set up species information
        ich = spc_info[name][0]
        smi = automol.inchi.smiles(ich)
        print("smiles: {}".format(smi), "inchi: {}".format(ich))

        # get and store the high level energies
        print('name in high level energy routine')
        print(name)
        for hl_idx, _ in enumerate(run_high_levels):
            min_ene = moldr.driver.get_high_level_energy(
                spc_info=spc_info[name],
                theory_low_level=ref_high_level,
                theory_high_level=run_high_levels[hl_idx],
                save_prefix=save_prefix,)
            ene_hl[(name, hl_idx)] = min_ene

    spc_save_fs = autofile.fs.species(save_prefix)
    for tors_model, vib_model in spc_models:
    # rot_model = RRHO or 1DHR or MDHR or TAU
        print('tors_model:', tors_model)
        print('vib_model:', vib_model)
        for har_level, tors_level, vpt2_level in pf_levels:
            # get the zero-point energy for each species
            print('harmonic_level:', har_level)
            print('tors_level:', tors_level)
            print('vpt2_level:', vpt2_level)
            print('Calculating zpe')
            spc_zpe = {}
            is_atom = {}
            zero_energy_str = {}
            for name in spc_names:
                ich = spc_info[name][0]
                smi = automol.inchi.smiles(ich)
                print("smiles: {}".format(smi), "inchi: {}".format(ich))
                spc_save_fs.leaf.create(spc_info[name])
                spc_save_path = spc_save_fs.leaf.path(spc_info[name])

                spc_zpe[name], is_atom[name] = moldr.driver.get_zero_point_energy(
                    spc_info[name],
                    tors_model, vib_model,
                    har_level, tors_level, vpt2_level,
                    script_str=pf_script_str,
                    elec_levels=[[0., 1]], sym_factor=1.,
                    save_prefix=spc_save_path)
                print(name, spc_zpe[name])
                zpe_str = '{0:<8.2f}\n'.format(spc_zpe[name])
                if is_atom[name]:
                   zero_energy_str[name] = 'End'
                else:
                   zero_energy_str[name] = ' ZeroEnergy[kcal/mol] ' + zpe_str
                   zero_energy_str[name] += 'End'
                print('zero_energy_str:', zero_energy_str)

            pf_inp_str = {}
            print('finished zpe')
            # get the partition function for each species
            
            chemkin_poly_strs = ['' for i in range(len(run_high_levels))]
            for name in spc_names:
                # set up species information
                ich = spc_info[name][0]
                smi = automol.inchi.smiles(ich)
                mul = spc_info[name][2]
                print("smiles: {}".format(smi), "inchi: {}".format(ich))
                spc_save_fs.leaf.create(spc_info[name])
                spc_save_path = spc_save_fs.leaf.path(spc_info[name])

                # to be generalized
                sym_factor = 1.
                elec_levels = [[0., mul]]
                if (ich, mul) in elc_deg_dct:
                    elec_levels = elc_deg_dct[(ich, mul)]

                # cycle through the low levels generating partition functions  for each
                spc_str = moldr.driver.species_block(
                    spc_info=spc_info[name],
                    tors_model=tors_model,
                    vib_model=vib_model,
                    har_level=har_level,
                    tors_level=tors_level,
                    vpt2_level=vpt2_level,
                    script_str=projrot_script_str,
                    elec_levels=elec_levels,
                    sym_factor=sym_factor,
                    save_prefix=spc_save_path,
                    )

                # create a messpf input file
                global_pf_str = mess_io.writer.write_global_pf(
                    [], temp_step, ntemps, rel_temp_inc=0.001,
                    atom_dist_min=0.6)
                print(global_pf_str)
                spc_head_str = 'Species ' + name
                print(spc_head_str)
                pf_inp_str[name] = '\n'.join(
                    [global_pf_str, spc_head_str,
                     spc_str, zero_energy_str[name], '\n'])
                print(spc_str)
                print(pf_inp_str[name])

                orb_restr = moldr.util.orbital_restriction(
                        spc_info[name], tors_level)
                tors_levelp = tors_level[1:3]
                tors_levelp.append(orb_restr)

                thy_save_fs = autofile.fs.theory(spc_save_path)
                thy_save_fs.leaf.create(tors_levelp)
                thy_save_path = thy_save_fs.leaf.path(tors_levelp)
                bld_locs = ['PF', 0]
                bld_save_fs = autofile.fs.build(thy_save_path)
                bld_save_fs.leaf.create(bld_locs)
                pf_path = bld_save_fs.leaf.path(bld_locs)
                print('Build Path for Partition Functions')
                print(pf_path)

                # run messpf
                with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
                    pf_file.write(pf_inp_str[name])
                moldr.util.run_script(pf_script_str, pf_path)
                print('finished partition function')

                formula = thermo.util.inchi_formula(ich)
                print('\nformula:')
                print(formula)

                # Get atom count dictionary
                atom_dict = thermo.util.get_atom_counts_dict(formula)
                print('\natom dict:')
                print(atom_dict)

                # Get the list of the basis
                spc_bas = thermo.heatform.get_reduced_basis(spc_ref_ich, formula)
                print('\nbasis:')
                print(spc_bas)

                # Get the coefficients for the balanced heat-of-formation eqn
                coeff = thermo.heatform.calc_coefficients(spc_bas, atom_dict)
                print('\ncoeff:')
                print(coeff)

                # prepare NASA polynomials
                nasa_inp_str = ('nasa')
                bld_locs = ['NASA_POLY', 0]
                bld_save_fs.leaf.create(bld_locs)
                nasa_path = bld_save_fs.leaf.path(bld_locs)
                print('NASA build path')
                print(nasa_path)

                h_basis = []
                for hl_idx, _ in enumerate(run_high_levels):
                    for  name_ref in spc_ref_names:
                        ich = spc_info[name_ref][0]
                        if ich in spc_bas:
                            tmp = ene_hl[(name_ref, hl_idx)] + spc_zpe[name_ref]/EH2KCAL
                            h_basis.append(tmp)
                    print('\ne_basis:')
                    print(h_basis)

                    # Get the 0 K heat of formation
                    spc_ene = ene_hl[(name, hl_idx)] + spc_zpe[name]/EH2KCAL
                    h0form = thermo.heatform.calc_hform_0k(spc_ene, h_basis, spc_bas, coeff, ref_set='ATcT')
                    print('h0form = ',h0form)

                    # need to change back to starting directory after running thermp and pac99 or rest of code is confused
                    starting_path = os.getcwd()
                    os.chdir(nasa_path)

                    # Write thermp input file
                    ENTHALPYT = 0.
                    BREAKT = 1000.
                    thermo.runner.write_thermp_input(
                        formula=formula,
                        deltaH=h0form,
                        enthalpyT=ENTHALPYT,
                        breakT=BREAKT,
                        thermp_file_name='thermp.dat')

                    # Run thermp
                    thermo.runner.run_thermp(
                        pf_path=pf_path,
                        thermp_path=nasa_path,
                        thermp_file_name='thermp.dat',
                        pf_file_name='pf.dat'
                        )

                    # Run pac99
                    print('formula test')
                    print(formula)
                    print(nasa_path)
                    thermo.runner.run_pac99(nasa_path, formula)

                    with open(os.path.join(nasa_path, 'thermp.out'), 'r') as thermp_outfile:
                        thermp_out_str = thermp_outfile.read()

                    with open(os.path.join(nasa_path, formula+'.c97'), 'r') as pac99_file:
                        pac99_str = pac99_file.read()

                    # Get the pac99 polynomial
                    pac99_poly_str = thermo.nasapoly.get_pac99_polynomial(pac99_str)
                    print('\nPAC99 Polynomial:')
                    print(pac99_poly_str)

                    # Convert the pac99 polynomial to chemkin polynomial
                    comment_str = '! tors model: {0}\n'.format(tors_model)
                    comment_str += '! vib model: {0}\n'.format(vib_model)
                    comment_str += '! har level: {0}{1}/{2}\n'.format(
                        har_level[3], har_level[1], har_level[2])
                    comment_str += '! tors level: {0}{1}/{2}\n'.format(
                        tors_level[3], tors_level[1], tors_level[2])
                    comment_str += '! vpt2 level: {0}{1}/{2}\n'.format(
                        vpt2_level[3], vpt2_level[1], vpt2_level[2])
                    comment_str += '! ref level for energy: {0}{1}/{2}\n'.format(
                        ref_high_level[3], ref_high_level[1], ref_high_level[2])
                    comment_str += '! energy level: {0}/{1}\n'.format(
                        run_high_levels[hl_idx][1], run_high_levels[hl_idx][2])

                    chemkin_poly_str = thermo.nasapoly.convert_pac_to_chemkin(
                        name, atom_dict, comment_str, pac99_poly_str)
                    chemkin_set_str = thermo.nasapoly.convert_pac_to_chemkin(
                        name, atom_dict, '', pac99_poly_str)
                    print('\nCHEMKIN Polynomial:')
                    print(chemkin_poly_str)
                    print(hl_idx)
                    if chemkin_poly_strs[hl_idx] == '':
                        chemkin_poly_strs[hl_idx] += chemkin_poly_str
                    else:
                        chemkin_poly_strs[hl_idx] += chemkin_set_str
                    print('startig_path in thermo')
                    print(starting_path)
                    ckin_path = ''.join([starting_path, '/ckin'])
                    print(ckin_path)
                    os.chdir(starting_path)
                    with open(os.path.join(nasa_path, formula+'.ckin'), 'w') as nasa_file:
                        nasa_file.write(chemkin_poly_str)
                    with open(os.path.join(ckin_path, formula+'.ckin'), 'w') as nasa_file:
                        nasa_file.write(chemkin_poly_str)

            for hl_idx, _ in enumerate(run_high_levels):
                hl_idx_str = str(hl_idx)
                with open(os.path.join(ckin_path, 'SPECIES'+hl_idx_str+'.ckin'), 'w') as nasa_file:
                    nasa_file.write(chemkin_poly_strs[hl_idx])


def reaction_rates(
        rct_names_lst,
        prd_names_lst,
        smi_dct,
        chg_dct,
        mul_dct,
        elc_deg_dct,
        temperatures,
        pressures,
        etsfr_par,
        run_opt_levels,
        ref_high_level,
        run_high_levels,
        spc_models,
        pf_levels,
        save_prefix,
        projrot_script_str,
        rate_script_str,
    ):
    """Evaluate the reaction rates for given ab initio levels and theory models
    """

    idx = 0
    exp_factor = etsfr_par[idx]
    idx += 1

    exp_power = etsfr_par[idx]
    idx += 1

    exp_cutoff = etsfr_par[idx]
    idx += 1

    eps1 = etsfr_par[idx]
    idx += 1

    eps2 = etsfr_par[idx]
    idx += 1

    sig1 = etsfr_par[idx]
    idx += 1

    sig2 = etsfr_par[idx]
    idx += 1

    mass1 = etsfr_par[idx]
    idx += 1

    mass2 = etsfr_par[idx]
    idx += 1

    energy_trans_str = mess_io.writer.write_energy_transfer(
        exp_factor, exp_power, exp_cutoff, eps1, eps2, sig1, sig2, mass1, mass2)
    print(energy_trans_str)

    for rct_names, prd_names in zip(rct_names_lst, prd_names_lst):
        # print the CHEMKIN reaction name for reference
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
        print()
        print('Mess Input for')
        print("Reaction: {}".format(rxn_name))

        # determine inchis, charges, and multiplicities

        rct_smis = list(map(smi_dct.__getitem__, rct_names))
        prd_smis = list(map(smi_dct.__getitem__, prd_names))
        rct_ichs = list(map(automol.smiles.inchi, rct_smis))
        prd_ichs = list(map(automol.smiles.inchi, prd_smis))
        rct_chgs = list(map(chg_dct.__getitem__, rct_names))
        prd_chgs = list(map(chg_dct.__getitem__, prd_names))
        rct_muls = list(map(mul_dct.__getitem__, rct_names))
        prd_muls = list(map(mul_dct.__getitem__, prd_names))
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])

        ts_mul = automol.mult.ts.low(rct_muls, prd_muls)
        ts_chg = sum(rct_chgs)
        print('ts_chg test:',ts_chg)
        ts_info = ('', ts_chg, ts_mul)

# header section
        header_str = mess_io.writer.write_global_reaction(temperatures, pressures)
        print(header_str)

# energy transfer section

    for tors_model, vib_model in spc_models:
    # rot_model = RRHO or 1DHR or MDHR or TAU
        print('tors_model:',tors_model)
        print('vib_model:',vib_model)
        for har_level, tors_level, vpt2_level in pf_levels:
            # get the zero-point energy for each species
            print('harmonic_level:',har_level)
            print('tors_level:',tors_level)
            print('vpt2_level:',vpt2_level)
        for opt_level_idx, _ in enumerate(run_opt_levels):
#        for prog, method, basis in run_opt_levels:
            prog = run_opt_levels[opt_level_idx][0]
            SCRIPT_STR, OPT_SCRIPT_STR, KWARGS, OPT_KWARGS = (
                moldr.util.run_qchem_par(prog))

            ts_orb_restr = moldr.util.orbital_restriction(
                ts_info, run_opt_levels[opt_level_idx])
            thy_level = run_opt_levels[opt_level_idx][1:3]
            thy_level.append(ts_orb_restr)

            # check direction of reaction
            rxn_ichs = [rct_ichs, prd_ichs]
            rxn_chgs = [rct_chgs, prd_chgs]
            rxn_muls = [rct_muls, prd_muls]
            rxn_exo = moldr.util.reaction_energy(
                save_prefix, rxn_ichs, rxn_chgs, rxn_muls, run_opt_levels[opt_level_idx])
            print(rxn_exo)
            if rxn_exo > 0:
                rct_ichs, prd_ichs = prd_ichs, rct_ichs
                rct_chgs, prd_chgs = prd_chgs, rct_chgs
                rct_muls, prd_muls = prd_muls, rct_muls
                print('ts search will be performed in reverse direction')

            # construct the filesystem
            rxn_ichs = [rct_ichs, prd_ichs]
            rxn_chgs = [rct_chgs, prd_chgs]
            rxn_muls = [rct_muls, prd_muls]

            # set up the filesystem
            is_rev = autofile.system.reaction_is_reversed(
                rxn_ichs, rxn_chgs, rxn_muls)
            rxn_ichs, rxn_chgs, rxn_muls = autofile.system.sort_together(
                rxn_ichs, rxn_chgs, rxn_muls)
            print(" - The reaction direction is {}"
                  .format('backward' if is_rev else 'forward'))

            ts_mul = automol.mult.ts.low(rct_muls, prd_muls)
            ts_chg = 0
            for rct_chg in rct_chgs:
                ts_chg += rct_chg
            ts_info = ['', ts_chg, ts_mul]

            rxn_ichs = tuple(map(tuple, rxn_ichs))
            rxn_chgs = tuple(map(tuple, rxn_chgs))
            rxn_muls = tuple(map(tuple, rxn_muls))
            print('rxn_save test', rxn_ichs, rxn_chgs, rxn_muls, ts_mul)
            print(save_prefix)
            rxn_save_fs = autofile.fs.reaction(save_prefix)
#            rxn_save_fs.leaf.create([rxn_ichs, rxn_chgs, rxn_muls, ts_mul])
            rxn_save_path = rxn_save_fs.leaf.path(
                [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

            orb_restr = moldr.util.orbital_restriction(
                ts_info, run_opt_levels[opt_level_idx])
            ref_level = run_opt_levels[opt_level_idx][1:3]
            ref_level.append(orb_restr)
            print('ref level test:', ref_level)
            thy_save_fs = autofile.fs.theory(rxn_save_path)
            thy_save_fs.leaf.create(ref_level)
            thy_save_path = thy_save_fs.leaf.path(ref_level)
            spc_save_fs = autofile.fs.species(save_prefix)

            # cycle over reactant and product species
            # check if unimolecular or bimolecular species
            print('Unimolecular or bimolecular test')
            indxw = 0
            indxp = 0
            well_str = ''
            bim_str = ''
            # first reactants
            if not rct_ichs[1]:
                # well
                print(rct_ichs)
                indxw += 1
                spc_info = (rct_ichs[0], rct_chgs[0], rct_muls[0])
                spc_save_fs.leaf.create(spc_info)
                spc_save_path = spc_save_fs.leaf.path(spc_info)
                well_data = moldr.driver.species_block(
                    spc_info=spc_info,
                    tors_model=tors_model,
                    vib_model=vib_model,
                    har_level=har_level,
                    tors_level=tors_level,
                    vpt2_level=vpt2_level,
                    script_str=projrot_script_str,
                    elec_levels=[[0., 1]], sym_factor=1.,
                    save_prefix=spc_save_path,
                    )
                well_label = 'W'+str(indxw)
                reac_label = well_label
                zero_energy = 0.0
                well_str += mess_io.writer.write_well(well_label, well_data, zero_energy)
                print('well_str:', well_str)
                well_data.replace()
            else:
                # bimolecular
                indxp += 1
                spc_data = ['', '']
                spc_label = ['', '']
                bimol_label = 'P'+str(indxp)
                reac_label = bimol_label
                sym_factor = 1.
                for idx, (rct_ich, rct_chg, rct_mul) in enumerate(zip(rct_ichs, rct_chgs, rct_muls)):
                    spc_label[idx] = automol.inchi.smiles(rct_ich)
                    spc_info = (rct_ich, rct_chg, rct_mul)
                    spc_save_fs.leaf.create(spc_info)
                    spc_save_path = spc_save_fs.leaf.path(spc_info)
                    spc_data[idx] = moldr.driver.species_block(
                        spc_info=spc_info,
                        tors_model=tors_model,
                        vib_model=vib_model,
                        har_level=har_level,
                        tors_level=tors_level,
                        vpt2_level=vpt2_level,
                        script_str=projrot_script_str,
                        elec_levels=[[0., 1]], sym_factor=1.,
                        save_prefix=spc_save_path,
                        )
                ground_energy = 0.0
                bim_str += mess_io.writer.write_bimolecular(
                    bimol_label, spc_label[0], spc_data[0],
                    spc_label[1], spc_data[1], ground_energy)
                print('bim str:', bim_str)

            # now products
            if not prd_ichs[1]:
                # well
                print(prd_ichs)
                indxw += 1
                for prd_ich, prd_chg, prd_mul in zip(prd_ichs, prd_chgs, prd_muls):
                    spc_info = (prd_ich, prd_chg, prd_mul)
                    spc_save_fs.leaf.create(spc_info)
                    spc_save_path = spc_save_fs.leaf.path(spc_info)
                    well_data = moldr.driver.species_block(
                        spc_info=spc_info,
                        tors_model=tors_model,
                        vib_model=vib_model,
                        har_level=har_level,
                        tors_level=tors_level,
                        vpt2_level=vpt2_level,
                        script_str=projrot_script_str,
                        elec_levels=[[0., 1]], sym_factor=1.,
                        save_prefix=spc_save_path,
                        )
                well_label = 'W'+str(indxw)
                prod_label = well_label
                zero_energy = 0.0
                well_str += mess_io.writer.write_well(well_label, well_data, zero_energy)
                print('well_str:', well_str)
    # write W_indxw 
            else:
                # bimolecular
                indxp += 1
                spc_data = ['', '']
                spc_label = ['', '']
                bimol_label = 'P'+str(indxp)
                prod_label = bimol_label
                for idx, (prd_ich, prd_chg, prd_mul) in enumerate(zip(prd_ichs, prd_chgs, prd_muls)):
                    spc_info = (prd_ich, prd_chg, prd_mul)
                    spc_save_fs.leaf.create(spc_info)
                    spc_save_path = spc_save_fs.leaf.path(spc_info)
                    spc_data[idx] = moldr.driver.species_block(
                        spc_info=spc_info,
                        tors_model=tors_model,
                        vib_model=vib_model,
                        har_level=har_level,
                        tors_level=tors_level,
                        vpt2_level=vpt2_level,
                        script_str=projrot_script_str,
                        elec_levels=[[0., 1]], sym_factor=1.,
                        save_prefix=spc_save_path,
                        )
                    spc_label[idx] = automol.inchi.smiles(prd_ich)
                ground_energy = 0.0
                bim_str += mess_io.writer.write_bimolecular(
                    bimol_label, spc_label[0], spc_data[0],
                    spc_label[1], spc_data[1], ground_energy)
                print('bim str:', bim_str)

            elec_levels = [[0., ts_mul]]
            symfactor = 1.
#            prd_info=[prd_ichs[i], prd_chgs[i], prd_muls[i]]
            # is it OK to use product info??
#            spc_save_path = moldr.util.species_path(
#                prd_ichs[i], prd_chgs[i], prd_muls[i], save_prefix)
#            thy_save_path = moldr.util.theory_path(method, basis, ts_orb_restr, spc_save_path)
            print('thy_save_path:', thy_save_path)
            ts_data_str = moldr.driver.species_block(
                spc_info=ts_info,
                tors_model=tors_model,
                vib_model=vib_model,
                har_level=har_level,
                tors_level=tors_level,
                vpt2_level=vpt2_level,
                script_str=projrot_script_str,
                elec_levels=[[0., 1]], sym_factor=1.,
                save_prefix=thy_save_path,
                )
            ts_label = 'B1'

            ts_str = mess_io.writer.write_ts_sadpt(ts_label, reac_label, prod_label, ts_data_str, zero_energy)
            print(ts_str)

#                    ts_sadpt_writer
            core = mess_io.writer.write_core_rigidrotor(geo, symfactor)
            molecule_section_str1 = mess_io.writer.write_molecule(
                core, freqs, zpe, elec_levels, hind_rot='',
            )
            print(molecule_section_str1)
    #
            mess_inp_str = '/n'.join([header_str, energy_trans_str, well_str, bim_str, ts_str])

            bld_locs = ['MESS', 0]
            bld_save_fs = autofile.fs.build(thy_save_path)
            bld_save_fs.leaf.create(bld_locs)
            mess_path = bld_save_fs.leaf.path(bld_locs)
            print('Build Path for MESS rate files:')
            print(mess_path)
            with open(os.path.join(mess_path, 'mess.inp'), 'w') as mess_file:
                mess_file.write(mess_inp_str[name])
            moldr.util.run_script(rate_script_str, mess_path)
                
#            rct_block_str = []
#            for i, _ in enumerate(rct_ichs):
#                elec_levels = [[0., mul]]
#                if (ich, mul) in elc_deg_dct:
#                    elec_levels = elc_deg_dct[(ich, mul)]
#                symfactor = 1.
#                rct_info=[rct_ichs[i], rct_chgs[i], rct_muls[i]]
#                orb_restr = moldr.util.orbital_restriction(rct_muls[i], RESTRICT_OPEN_SHELL)
#                spc_save_path = moldr.util.species_path(
#                    rct_ichs[i], rct_chgs[i], rct_muls[i], save_prefix)
#                thy_save_path = moldr.util.theory_path(method, basis, orb_restr, spc_save_path)
#                rct_block_str.append(moldr.driver.species_block(
#                        spc_info=rct_info,
#                        tors_model=tors_model,
#                        vib_model=vib_model,
#                        har_level=har_level,
#                        tors_level=tors_level,
#                        vpt2_level=vpt2_level,
#                        orb_restr=orb_restr,
#                        script_str=projrot_script_str,
#                        elec_levels=elec_levels,
#                        sym_factor=sym_factor,
#                        save_prefix=thy_save_path,
#                        ))
#            indx_well = 0
#            indx_bim = 0
#            if rct_ichs[1]:
#                print(rct_ichs)
#                indx_well += 1
#                reac_label = 'W{}'.format(indx_well)
#
## write W_indxw
#            else:
#                indx_bim += 1
#
#                prd_block_str = []
#                for i, _ in enumerate(prd_ichs):
#                    elec_levels = [[0., mul]]
#                    if (ich, mul) in elc_deg_dct:
#                        elec_levels = elc_deg_dct[(ich, mul)]
#                    symfactor = 1.
#                    prd_info=[prd_ichs[i], prd_chgs[i], prd_muls[i]]
#                    orb_restr = moldr.util.orbital_restriction(prd_muls[i], RESTRICT_OPEN_SHELL)
#                    spc_save_path = moldr.util.species_path(
#                            prd_ichs[i], prd_chgs[i], prd_muls[i], save_prefix)
#                    thy_save_path = moldr.util.theory_path(method, basis, orb_restr, spc_save_path)
#                    prd_block_str.append(moldr.driver.species_block(
#                            spc_info=prd_info,
#                            tors_model=tors_model,
#                            vib_model=vib_model,
#                            har_level=har_level,
#                            tors_level=tors_level,
#                            vpt2_level=vpt2_level,
#                            script_str=projrot_script_str,
#                            elec_levels=elec_levels,
#                            sym_factor=sym_factor,
#                            save_prefix=thy_save_path,
#                            ))


        # write P_indxp
        # cycle over vdw Species
        # cycle over transition states
        # for now - have only one TS

        #        if reaction_typ = addition:
        # reactants
        #            print(spc_str(rct1))
        #            print(spc_str(rct2))
        # well
        #            print(spc_str(prod1))
        #        if reaction_typ = abstraction:
        # reactants
        #            print(spc_str(rct1))
        #            print(spc_str(rct2))
        # vdw
        #           for vdw_species in ...
        #                   print(spc_str(vdwi))
        # products
        #            print(spc_str(prod1))
        #            print(spc_str(prod2))

        #        if reaction_typ = abstraction
        # reactants
        #            print(spc_str(rct1))
        #            print(spc_str(rct2))
        # vdw
        #            print(spc_str(vdw1))
        #            print(spc_str(vdw2))
        # products
        #            print(spc_str(prod1))
        #            print(spc_str(prod2))
        # ts

