""" reaction list test
"""
import os
from qcelemental import constants as qcc
import thermo
import automol
import autofile
import moldr
import mess_io.writer

ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')
PF_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                 "messpf pf.inp build.out >> stdout.log &> stderr.log")
PROJROT_SCRIPT_STR = ("#!/usr/bin/env bash\n"
                      "RPHt.exe")
ELC_DEG_DCT = {
    ('InChI=1S/B', 2): [[0., 2], [16., 4]],
    ('InChI=1S/C', 3): [[0., 1], [16.4, 3], [43.5, 5]],
    ('InChI=1S/N', 2): [[0., 6], [8., 4]],
    ('InChI=1S/O', 3): [[0., 5], [158.5, 3], [226.5, 1]],
    ('InChI=1S/F', 2): [[0., 4], [404.1, 2]],
    ('InChI=1S/Cl', 2): [[0., 4], [883.4, 2]],
    ('InChI=1S/Br', 2): [[0., 4], [685.2, 2]],
    ('InChI=1S/HO/h1H', 2): [[0., 2], [138.9, 2]],
    ('InChI=1S/NO/c1-2', 2): [[0., 2], [123.1, 2]],
    ('InChI=1S/O2/c1-2', 1): [[0., 2]]
}

def get_electronic_energy(spc_info, geo_theory, sp_theory, save_prefix):
    """ return the electronic energy for a specific species for a given level of theory
    """
    ene = moldr.pf.get_high_level_energy(
        spc_info=spc_info,
        thy_low_level=geo_theory,
        thy_high_level=sp_theory,
        save_prefix=save_prefix)
    return ene

def spc_energy(spc_ene, spc_zpe):
    """ return the sum of the electronic and zero point energies
    """
    spc_ene = spc_ene + spc_zpe/EH2KCAL
    return spc_ene

def basis_energy(spc_bas, spcdct):
    """ return the electronic + zero point energies for a set of species
    """
    h_basis = []
    for ich in spc_bas:
        for key in spcdct:
            if ich == spcdct[key]['ich']:
                tmp = spc_energy(spcdct[key]['ene'], spcdct[key]['zpe'])
                h_basis.append(tmp)
                break
    return h_basis

def get_coeff(spc, spcdct, spc_bas):
    """ return the coefficients for the expansion of a species in terms of a
    set of basis species
    """
    ich = spcdct[spc]['ich']
    formula = thermo.util.inchi_formula(ich)
    print('\nformula:')
    print(formula)

    # Get atom count dictionary
    atom_dict = thermo.util.get_atom_counts_dict(formula)
    print('\natom dict:')
    print(atom_dict)

    if len(spc_bas) == 1:
        if spcdct[spc]['ich'] == spc_bas[0]:
            coeff = [1]
            print('\ncoeff:')
            print(coeff)
    else:
        # Get the coefficients for the balanced heat-of-formation eqn
        coeff = thermo.heatform.calc_coefficients(spc_bas, atom_dict)
        print('\ncoeff:')
        print(coeff)
    return coeff

def get_hf0k(spc, spcdct, spc_bas):
    """ determine the 0 K heat of formation from the
    species dictionary and a set of references species
    """
    spc_ene = spc_energy(spcdct[spc]['ene'], spcdct[spc]['zpe'])
    h_basis = basis_energy(spc_bas, spcdct)
    coeff = get_coeff(spc, spcdct, spc_bas)

    # Get the 0 K heat of formation
    # ref_set should be a parameter for this routine
    print('spc_ene:', spc_ene)
    print('h_basis:', h_basis)
    print('coeff:', coeff)
    h0form = thermo.heatform.calc_hform_0k(spc_ene, h_basis, spc_bas, coeff, ref_set='ATcT')
    print('Hf0K test:', spc, h0form)
    return h0form

def get_zpe(spcdct, spc_info, spc_save_path, pf_levels, spc_model):
    """ return the zpe for a given species according a specified set of
    partition function levels
    """
    har_level, tors_level, vpt2_level = pf_levels
    tors_model, vib_model = spc_model
    # get the zero-point energy for each species
    # print('harmonic_level:', har_level)
    # print('tors_level:', tors_level)
    # print('vpt2_level:', vpt2_level)
    # print('Calculating zpe')
    spc_zpe = {}
    is_atom = {}
    zero_energy_str = {}
    ich = spcdct['ich']
    smi = automol.inchi.smiles(ich)
    print("smiles: {}".format(smi), "inchi: {}".format(ich))

    spc_zpe, is_atom = moldr.pf.get_zero_point_energy(
        spc_info,
        tors_model, vib_model,
        har_level, tors_level, vpt2_level,
        script_str=PF_SCRIPT_STR,
        elec_levels=[[0., 1]], sym_factor=1.,
        save_prefix=spc_save_path)
    zpe_str = '{0:<8.2f}\n'.format(spc_zpe)
    print(zpe_str)
    if is_atom:
        zero_energy_str = 'End'
    else:
        zero_energy_str = ' ZeroEnergy[kcal/mol] ' + zpe_str
        zero_energy_str += 'End'
    print('zero_energy_str:', zero_energy_str)

    print('finished zpe')
    return spc_zpe, zero_energy_str

def get_spcinput(spcdct, spc_info, spc_save_path, pf_levels, spc_model):
    """ set up the input string for a given species section in mess input
    """

    har_level, tors_level, vpt2_level = pf_levels
    tors_model, vib_model = spc_model

    # set up species information
    ich = spc_info[0]
    smi = automol.inchi.smiles(ich)
    mul = spc_info[2]
    print("smiles: {}".format(smi), "inchi: {}".format(ich))

    sym_factor = 1.
    elec_levels = [[0., mul]]
    if 'sym' in spcdct:
        sym_factor = spcdct['sym']
    if 'elec_levs' in spcdct:
        elec_levels = spcdct['elec_levs']
    if (ich, mul) in  ELC_DEG_DCT:
        elec_levels = ELC_DEG_DCT[(ich, mul)]

    # cycle through the low levels generating partition functions  for each
    spc_str = moldr.pf.species_block(
        spc_info=spc_info,
        tors_model=tors_model,
        vib_model=vib_model,
        har_level=har_level,
        tors_level=tors_level,
        vpt2_level=vpt2_level,
        script_str=PROJROT_SCRIPT_STR,
        elec_levels=elec_levels,
        sym_factor=sym_factor,
        save_prefix=spc_save_path,
        )
    print(spc_str)
    return spc_str

def get_pfheader(temp_step, ntemps):
    """ prepare partition function header string
    """
    global_pf_str = mess_io.writer.global_pf(
        [], temp_step, ntemps, rel_temp_inc=0.001,
        atom_dist_min=0.6)
    print(global_pf_str)
    return global_pf_str

def get_thermo_paths(spc_save_path, spc_info, har_level):
    """ set up the path for saving the pf input and output
    currently using the harmonic theory directory for this because
    there is no obvious place to save this information for a random
    assortment of har_level, tors_level, vpt2_level
    """
    orb_restr = moldr.util.orbital_restriction(
        spc_info, har_level)
    har_levelp = har_level[1:3]
    har_levelp.append(orb_restr)

    thy_save_fs = autofile.fs.theory(spc_save_path)
    thy_save_fs.leaf.create(har_levelp)
    thy_save_path = thy_save_fs.leaf.path(har_levelp)
    bld_locs = ['PF', 0]
    bld_save_fs = autofile.fs.build(thy_save_path)
    bld_save_fs.leaf.create(bld_locs)
    pf_path = bld_save_fs.leaf.path(bld_locs)

    # prepare NASA polynomials
    bld_locs = ['NASA_POLY', 0]
    bld_save_fs.leaf.create(bld_locs)
    nasa_path = bld_save_fs.leaf.path(bld_locs)

    print('Build Path for Partition Functions')
    print(pf_path)
    print('NASA build path')
    print(nasa_path)
    return pf_path, nasa_path

def get_pfinput(spc, spc_str, global_pf_str, zpe_str):
    """ prepare the full pf input string for running messpf
    """

    # create a messpf input file
    spc_head_str = 'Species ' + spc
    pf_inp_str = '\n'.join(
        [global_pf_str, spc_head_str,
         spc_str, zpe_str, '\n'])
    print(spc_str)
    print(pf_inp_str)
    return pf_inp_str

def write_pfinput(pf_inp_str, pf_path):
    """ write the pf.inp file
    """
    # run messpf
    with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
        pf_file.write(pf_inp_str)

def run_pf(pf_path, pf_script_str=PF_SCRIPT_STR):
    """ run messpf
    """
    moldr.util.run_script(pf_script_str, pf_path)
    print('finished partition function')

def go_to_path(path):
    """ change directory to path and return the original working directory
    """
    starting_path = os.getcwd()
    os.chdir(path)
    return starting_path

def return_to_path(path):
    """ change directory to starting path
    """
    os.chdir(path)

def prepare_path(path, loc):
    """ change directory to starting path, return chemkin path
    """
    new_path = os.path.join(path, loc)
    return new_path

def write_thermp_inp(spc_spcdct):
    """ write the thermp input file
    """
    ich = spc_spcdct['ich']
    h0form = spc_spcdct['Hf0K']
    formula = thermo.util.inchi_formula(ich)
    # Write thermp input file
    enthalpyt = 0.
    breakt = 1000.
    thermo.runner.write_thermp_input(
        formula=formula,
        delta_h=h0form,
        enthalpy_temp=enthalpyt,
        break_temp=breakt,
        thermp_file_name='thermp.dat')

def run_thermp(pf_path, nasa_path):
    """ run thermp to convert partition functions to thermochemical data
    """
    # Run thermp
    thermo.runner.run_thermp(
        pf_path=pf_path,
        thermp_path=nasa_path,
        thermp_file_name='thermp.dat',
        pf_file_name='pf.dat'
        )

def run_pac(spc, spc_spcdct, nasa_path, pf_levels, spc_model):
    """ run pac99 to convert thermochemical data to nasa polynomials
    """
    ich = spc_spcdct['ich']
    formula = thermo.util.inchi_formula(ich)
    atom_dict = thermo.util.get_atom_counts_dict(formula)
    tors_model, vib_model = spc_model
    har_level, tors_level, vpt2_level, geo_level, sp_level = pf_levels

    # Run pac99
    print('formula test')
    print(formula)
    print(nasa_path)
    thermo.runner.run_pac99(nasa_path, formula)

    # this looks like it is not used - I wonder why it was there?
    # with open(os.path.join(nasa_path, 'thermp.out'), 'r') as thermp_outfile:
    #     thermp_out_str = thermp_outfile.read()

    with open(os.path.join(nasa_path, formula+'.c97'), 'r') as pac99_file:
        pac99_str = pac99_file.read()

    # Get the pac99 polynomial
    pac99_poly_str = thermo.nasapoly.get_pac99_polynomial(pac99_str)
    print('\nPAC99 Polynomial:')
    print(pac99_poly_str)

    # Convert the pac99 polynomial to chemkin polynomial
    if tors_model:
        comment_str = '! tors model: {0}\n'.format(tors_model)
    if vib_model:
        comment_str += '! vib model: {0}\n'.format(vib_model)
    if har_level:
        comment_str += '! har level: {0}{1}/{2}\n'.format(
            har_level[3], har_level[1], har_level[2])
    if tors_level:
        comment_str += '! tors level: {0}{1}/{2}\n'.format(
            tors_level[3], tors_level[1], tors_level[2])
    if vpt2_level:
        comment_str += '! vpt2 level: {0}{1}/{2}\n'.format(
            vpt2_level[3], vpt2_level[1], vpt2_level[2])
    if geo_level:
        comment_str += '! ref level for energy: {0}{1}/{2}\n'.format(
            geo_level[3], geo_level[1], geo_level[2])
    if sp_level:
        comment_str += '! energy level: {0}/{1}\n'.format(
            sp_level[1], sp_level[2])

    chemkin_poly_str = thermo.nasapoly.convert_pac_to_chemkin(
        spc, atom_dict, comment_str, pac99_poly_str)
    chemkin_set_str = thermo.nasapoly.convert_pac_to_chemkin(
        spc, atom_dict, '', pac99_poly_str)
    print('\nCHEMKIN Polynomial:')
    print(chemkin_poly_str)
    return chemkin_poly_str

def write_nasa_file(spc_spcdct, ckin_path, nasa_path, chemkin_poly_str):
    """ write out the nasa polynomials
    """
    ich = spc_spcdct['ich']
    formula = thermo.util.inchi_formula(ich)
    with open(os.path.join(nasa_path, formula+'.ckin'), 'w') as nasa_file:
        nasa_file.write(chemkin_poly_str)
    with open(os.path.join(ckin_path, formula+'.ckin'), 'w') as nasa_file:
        nasa_file.write(chemkin_poly_str)
    return
# create a messpf input file
#
#def species_thermo(
#        spc_names,
#        spc_info,
#        spc_ref_names,
#        elc_deg_dct,
#        temp_par,
#        ref_high_level,
#        run_high_levels,
#        spc_models,
#        pf_levels,
#        save_prefix,
#        projrot_script_str,
#        pf_script_str,
#        ):
#    """Evaluate the thermodynamics properties for all the species
#    """
#    # initializations
#    temp_step = temp_par[0]
#    ntemps = temp_par[1]
#    chemkin_poly_strs = {}
#    for hl_idx, _ in enumerate(run_high_levels):
#        chemkin_poly_strs[hl_idx] = ''
#    spc_ref_ich = []
#    for ref_name in spc_ref_names:
#        spc_ref_ich.append(spc_info[ref_name][0])
#
#    # generate and store the energies so that the reference energies can be used as needed
#    ene_hl = {}
#    print('Evaluating thermo for following species')
#    print(spc_names)
#    for name in spc_names:
#        # set up species information
#        ich = spc_info[name][0]
#        smi = automol.inchi.smiles(ich)
#        print("smiles: {}".format(smi), "inchi: {}".format(ich))
#
#        # get and store the high level energies
#        print('name in high level energy routine')
#        print(name)
#        for hl_idx, _ in enumerate(run_high_levels):
#
#    spc_save_fs = autofile.fs.species(save_prefix)
#    for tors_model, vib_model in spc_models:
#    # rot_model = RRHO or 1DHR or MDHR or TAU
#        print('tors_model:', tors_model)
#        print('vib_model:', vib_model)
#        for har_level, tors_level, vpt2_level in pf_levels:
#            # get the zero-point energy for each species
#            print('harmonic_level:', har_level)
#            print('tors_level:', tors_level)
#            print('vpt2_level:', vpt2_level)
#            print('Calculating zpe')
#            spc_zpe = {}
#            is_atom = {}
#            zero_energy_str = {}
#            for name in spc_names:
#                ich = spc_info[name][0]
#                smi = automol.inchi.smiles(ich)
#                print("smiles: {}".format(smi), "inchi: {}".format(ich))
#                spc_save_fs.leaf.create(spc_info[name])
#                spc_save_path = spc_save_fs.leaf.path(spc_info[name])
#
#                spc_zpe[name], is_atom[name] = moldr.pf.get_zero_point_energy(
#                    spc_info[name],
#                    tors_model, vib_model,
#                    har_level, tors_level, vpt2_level,
#                    script_str=pf_script_str,
#                    elec_levels=[[0., 1]], sym_factor=1.,
#                    save_prefix=spc_save_path)
#                print(name, spc_zpe[name])
#                zpe_str = '{0:<8.2f}\n'.format(spc_zpe[name])
#                if is_atom[name]:
#                   zero_energy_str[name] = 'End'
#                else:
#                   zero_energy_str[name] = ' ZeroEnergy[kcal/mol] ' + zpe_str
#                   zero_energy_str[name] += 'End'
#                print('zero_energy_str:', zero_energy_str)
#
#            pf_inp_str = {}
#            print('finished zpe')
#            # get the partition function for each species
#            
#            chemkin_poly_strs = ['' for i in range(len(run_high_levels))]
#            for name in spc_names:
#                # set up species information
#                ich = spc_info[name][0]
#                smi = automol.inchi.smiles(ich)
#                mul = spc_info[name][2]
#                print("smiles: {}".format(smi), "inchi: {}".format(ich))
#                spc_save_fs.leaf.create(spc_info[name])
#                spc_save_path = spc_save_fs.leaf.path(spc_info[name])
#
#                # to be generalized
#                sym_factor = 1.
#                elec_levels = [[0., mul]]
#                if (ich, mul) in elc_deg_dct:
#                    elec_levels = elc_deg_dct[(ich, mul)]
#
#                # cycle through the low levels generating partition functions  for each
#                spc_str = moldr.pf.species_block(
#                    spc_info=spc_info[name],
#                    tors_model=tors_model,
#                    vib_model=vib_model,
#                    har_level=har_level,
#                    tors_level=tors_level,
#                    vpt2_level=vpt2_level,
#                    script_str=projrot_script_str,
#                    elec_levels=elec_levels,
#                    sym_factor=sym_factor,
#                    save_prefix=spc_save_path,
#                    )
#
#                # create a messpf input file
#                global_pf_str = mess_io.writer.write_global_pf(
#                    [], temp_step, ntemps, rel_temp_inc=0.001,
#                    atom_dist_min=0.6)
#                print(global_pf_str)
#                spc_head_str = 'Species ' + name
#                print(spc_head_str)
#                pf_inp_str[name] = '\n'.join(
#                    [global_pf_str, spc_head_str,
#                     spc_str, zero_energy_str[name], '\n'])
#                print(spc_str)
#                print(pf_inp_str[name])
#
#                orb_restr = moldr.util.orbital_restriction(
#                        spc_info[name], tors_level)
#                tors_levelp = tors_level[1:3]
#                tors_levelp.append(orb_restr)
#
#                thy_save_fs = autofile.fs.theory(spc_save_path)
#                thy_save_fs.leaf.create(tors_levelp)
#                thy_save_path = thy_save_fs.leaf.path(tors_levelp)
#                bld_locs = ['PF', 0]
#                bld_save_fs = autofile.fs.build(thy_save_path)
#                bld_save_fs.leaf.create(bld_locs)
#                pf_path = bld_save_fs.leaf.path(bld_locs)
#                print('Build Path for Partition Functions')
#                print(pf_path)
#
#                # run messpf
#                with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
#                    pf_file.write(pf_inp_str[name])
#                moldr.util.run_script(pf_script_str, pf_path)
#                print('finished partition function')
#
#                formula = thermo.util.inchi_formula(ich)
#                print('\nformula:')
#                print(formula)
#
#                # Get atom count dictionary
#                atom_dict = thermo.util.get_atom_counts_dict(formula)
#                print('\natom dict:')
#                print(atom_dict)
#
#                # Get the list of the basis
#                spc_bas = thermo.heatform.get_reduced_basis(spc_ref_ich, formula)
#                print('\nbasis:')
#                print(spc_bas)
#
#                # Get the coefficients for the balanced heat-of-formation eqn
#                coeff = thermo.heatform.calc_coefficients(spc_bas, atom_dict)
#                print('\ncoeff:')
#                print(coeff)
#
#                # prepare NASA polynomials
#                nasa_inp_str = ('nasa')
#                bld_locs = ['NASA_POLY', 0]
#                bld_save_fs.leaf.create(bld_locs)
#                nasa_path = bld_save_fs.leaf.path(bld_locs)
#                print('NASA build path')
#                print(nasa_path)
#
#                h_basis = []
#                for hl_idx, _ in enumerate(run_high_levels):
#                    for  name_ref in spc_ref_names:
#                        ich = spc_info[name_ref][0]
#                        if ich in spc_bas:
#                            tmp = ene_hl[(name_ref, hl_idx)] + spc_zpe[name_ref]/EH2KCAL
#                            h_basis.append(tmp)
#                    print('\ne_basis:')
#                    print(h_basis)
#
#                    # Get the 0 K heat of formation
#                    spc_ene = ene_hl[(name, hl_idx)] + spc_zpe[name]/EH2KCAL
#                    h0form = thermo.heatform.calc_hform_0k(spc_ene, h_basis, spc_bas, coeff, ref_set='ATcT')
#                    print('h0form = ',h0form)
#
#                    # need to change back to starting directory after running thermp and pac99 or rest of code is confused
#                    starting_path = os.getcwd()
#                    os.chdir(nasa_path)
#
#                    # Write thermp input file
#                    ENTHALPYT = 0.
#                    BREAKT = 1000.
#                    thermo.runner.write_thermp_input(
#                        formula=formula,
#                        deltaH=h0form,
#                        enthalpyT=ENTHALPYT,
#                        breakT=BREAKT,
#                        thermp_file_name='thermp.dat')
#
#                    # Run thermp
#                    thermo.runner.run_thermp(
#                        pf_path=pf_path,
#                        thermp_path=nasa_path,
#                        thermp_file_name='thermp.dat',
#                        pf_file_name='pf.dat'
#                        )
#
#                    # Run pac99
#                    print('formula test')
#                    print(formula)
#                    print(nasa_path)
#                    thermo.runner.run_pac99(nasa_path, formula)
#
#                    with open(os.path.join(nasa_path, 'thermp.out'), 'r') as thermp_outfile:
#                        thermp_out_str = thermp_outfile.read()
#
#                    with open(os.path.join(nasa_path, formula+'.c97'), 'r') as pac99_file:
#                        pac99_str = pac99_file.read()
#
#                    # Get the pac99 polynomial
#                    pac99_poly_str = thermo.nasapoly.get_pac99_polynomial(pac99_str)
#                    print('\nPAC99 Polynomial:')
#                    print(pac99_poly_str)
#
#                    # Convert the pac99 polynomial to chemkin polynomial
#                    comment_str = '! tors model: {0}\n'.format(tors_model)
#                    comment_str += '! vib model: {0}\n'.format(vib_model)
#                    comment_str += '! har level: {0}{1}/{2}\n'.format(
#                        har_level[3], har_level[1], har_level[2])
#                    comment_str += '! tors level: {0}{1}/{2}\n'.format(
#                        tors_level[3], tors_level[1], tors_level[2])
#                    comment_str += '! vpt2 level: {0}{1}/{2}\n'.format(
#                        vpt2_level[3], vpt2_level[1], vpt2_level[2])
#                    comment_str += '! ref level for energy: {0}{1}/{2}\n'.format(
#                        ref_high_level[3], ref_high_level[1], ref_high_level[2])
#                    comment_str += '! energy level: {0}/{1}\n'.format(
#                        run_high_levels[hl_idx][1], run_high_levels[hl_idx][2])
#
#                    chemkin_poly_str = thermo.nasapoly.convert_pac_to_chemkin(
#                        name, atom_dict, comment_str, pac99_poly_str)
#                    chemkin_set_str = thermo.nasapoly.convert_pac_to_chemkin(
#                        name, atom_dict, '', pac99_poly_str)
#                    print('\nCHEMKIN Polynomial:')
#                    print(chemkin_poly_str)
#                    print(hl_idx)
#                    if chemkin_poly_strs[hl_idx] == '':
#                        chemkin_poly_strs[hl_idx] += chemkin_poly_str
#                    else:
#                        chemkin_poly_strs[hl_idx] += chemkin_set_str
#                    print('startig_path in thermo')
#                    print(starting_path)
#                    ckin_path = ''.join([starting_path, '/ckin'])
#                    print(ckin_path)
#                    os.chdir(starting_path)
#                    with open(os.path.join(nasa_path, formula+'.ckin'), 'w') as nasa_file:
#                        nasa_file.write(chemkin_poly_str)
#                    with open(os.path.join(ckin_path, formula+'.ckin'), 'w') as nasa_file:
#                        nasa_file.write(chemkin_poly_str)
#
#            for hl_idx, _ in enumerate(run_high_levels):
#                hl_idx_str = str(hl_idx)
#                with open(os.path.join(ckin_path, 'SPECIES'+hl_idx_str+'.ckin'), 'w') as nasa_file:
#                    nasa_file.write(chemkin_poly_strs[hl_idx])


