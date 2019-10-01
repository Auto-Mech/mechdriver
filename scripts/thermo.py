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
                      "RPHt.exe >& /dev/null")
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


def get_electronic_energy(spc_info, geo_theory, sp_theory, save_prefix, saddle=False):
    """ return the electronic energy for a specific species for a given level of theory
    """
    ene = moldr.pf.get_high_level_energy(
        spc_info=spc_info,
        thy_low_level=geo_theory,
        thy_high_level=sp_theory,
        save_prefix=save_prefix,
        saddle=saddle)
    return ene


def spc_energy(spc_ene, spc_zpe):
    """ return the sum of the electronic and zero point energies
    """
    spc_ene = spc_ene + spc_zpe/EH2KCAL
    return spc_ene


def basis_energy(spc_bas, spc_dct):
    """ return the electronic + zero point energies for a set of species
    """
    h_basis = []
    for ich in spc_bas:
        for key in spc_dct:
            if ich == spc_dct[key]['ich']:
                tmp = spc_energy(spc_dct[key]['ene'], spc_dct[key]['zpe'])
                h_basis.append(tmp)
                break
    return h_basis


def get_coeff(spc, spc_dct, spc_bas):
    """ return the coefficients for the expansion of a species in terms of a
    set of basis species
    """
    ich = spc_dct[spc]['ich']
    formula = automol.inchi.formula(ich)

    # Get atom count dictionary
    atom_dict = thermo.util.get_atom_counts_dict(formula)

    if len(spc_bas) == 1 and spc_dct[spc]['ich'] == spc_bas[0]:
        coeff = [1]
    else:
        # Get the coefficients for the balanced heat-of-formation eqn
        coeff = thermo.heatform.calc_coefficients(spc_bas, atom_dict)
    return coeff


def get_hf0k(spc, spc_dct, spc_bas):
    """ determine the 0 K heat of formation from the
    species dictionary and a set of references species
    """
    spc_ene = spc_energy(spc_dct[spc]['ene'], spc_dct[spc]['zpe'])
    h_basis = basis_energy(spc_bas, spc_dct)
    coeff = get_coeff(spc, spc_dct, spc_bas)

    # Get the 0 K heat of formation
    # ref_set should be a parameter for this routine
    h0form = thermo.heatform.calc_hform_0k(spc_ene, h_basis, spc_bas, coeff, ref_set='ATcT')
    return h0form


def get_zpe(spc, spc_info, spc_save_path, pf_levels, spc_model):
    """ return the zpe for a given species according a specified set of
    partition function levels
    """
    # get the zero-point energy for each species
    # print('harmonic_level:', har_level)
    # print('tors_level:', tors_level)
    # print('vpt2_level:', vpt2_level)
    # print('Calculating zpe')
    spc_zpe = {}
    is_atom = {}
    zero_energy_str = {}

    spc_zpe, is_atom = moldr.pf.get_zero_point_energy(
        spc, spc_info, pf_levels, spc_model,
        pf_script_str=PF_SCRIPT_STR,
        elec_levels=[[0., 1]], sym_factor=1.,
        save_prefix=spc_save_path)
    zpe_str = '{0:<8.2f}\n'.format(spc_zpe)
    if is_atom:
        zero_energy_str = 'End'
    else:
        zero_energy_str = ' ZeroEnergy[kcal/mol] ' + zpe_str
        zero_energy_str += 'End'

    return spc_zpe, zero_energy_str


def get_spc_input(spc, spc_dct_i, spc_info, spc_save_path, pf_levels, spc_model):
    """ set up the input string for a given species section in mess input
    """

    # set up species information
    ich = spc_info[0]
    smi = automol.inchi.smiles(ich)
    print("smiles: {}".format(smi), "inchi: {}".format(ich))

    # generate the partition function
    spc_str = moldr.pf.species_block(
        spc=spc,
        spc_dct_i=spc_dct_i,
        spc_info=spc_info,
        spc_model=spc_model,
        pf_levels=pf_levels,
        projrot_script_str=PROJROT_SCRIPT_STR,
        save_prefix=spc_save_path,
        )
    print(spc_str)
    return spc_str


def get_pf_header(temp_step, ntemps):
    """ prepare partition function header string
    """
    global_pf_str = mess_io.writer.global_pf(
        [], temp_step, ntemps, rel_temp_inc=0.001,
        atom_dist_min=0.6)
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
    return pf_path, nasa_path


def get_pf_input(spc, spc_str, global_pf_str, zpe_str):
    """ prepare the full pf input string for running messpf
    """

    # create a messpf input file
    spc_head_str = 'Species ' + spc
    pf_inp_str = '\n'.join(
        [global_pf_str, spc_head_str,
         spc_str, zpe_str, '\n'])
    return pf_inp_str


def write_pf_input(pf_inp_str, pf_path):
    """ write the pf.inp file
    """
    # run messpf
    with open(os.path.join(pf_path, 'pf.inp'), 'w') as pf_file:
        pf_file.write(pf_inp_str)


def run_pf(pf_path, pf_script_str=PF_SCRIPT_STR):
    """ run messpf
    """
    moldr.util.run_script(pf_script_str, pf_path)


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


def write_thermp_inp(spc_dct_i):
    """ write the thermp input file
    """
    ich = spc_dct_i['ich']
    h0form = spc_dct_i['Hfs'][0]
    formula = automol.inchi.formula(ich)
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
    with open('thermp.out', 'r') as f:
        lines = f.readlines()
    line = lines[-1]
    hf298k = line.split()[-1]
    return hf298k


def run_pac(spc_dct_i, nasa_path):
    """ run pac99 to convert thermochemical data to nasa polynomials
    """
    ich = spc_dct_i['ich']
    formula = automol.inchi.formula(ich)

    # Run pac99
    thermo.runner.run_pac99(nasa_path, formula)

    # this looks like it is not used - I wonder why it was there?
    # with open(os.path.join(nasa_path, 'thermp.out'), 'r') as thermp_outfile:
    #     thermp_out_str = thermp_outfile.read()

    with open(os.path.join(nasa_path, formula+'.c97'), 'r') as pac99_file:
        pac99_str = pac99_file.read()

    # Get the pac99 polynomial
    pac99_poly_str = thermo.nasapoly.get_pac99_polynomial(pac99_str)
    #print('\nPAC99 Polynomial:')
    #print('pac99_poly_str test:', pac99_poly_str)

    return pac99_poly_str


def run_ckin_header(pf_info, ref_info, spc_model):
    """ prepare chemkin header info and convert pac 99 format to chemkin format
    """
    tors_model, vib_model, _ = spc_model
    har_info, tors_info, vpt2_info, _, sp_str, = pf_info
    har_ref_info, tors_ref_info, vpt2_ref_info, _ = ref_info

    # Convert the pac99 polynomial to chemkin polynomial

    chemkin_header_str = '! tors model: {0}\n'.format(tors_model)
    chemkin_header_str += '! vib model: {0}\n'.format(vib_model)
    if har_info:
        chemkin_header_str += '! har level: {}{}/{}//{}{}/{}\n'.format(
            har_info[3], har_info[1], har_info[2],
            har_ref_info[3], har_ref_info[1], har_ref_info[2])
    if tors_info:
        chemkin_header_str += '! tors level: {}{}/{}//{}{}/{}\n'.format(
            tors_info[3], tors_info[1], tors_info[2],
            tors_ref_info[3], tors_ref_info[1], tors_ref_info[2])
    if vpt2_info:
        chemkin_header_str += '! vpt2 level: {}{}/{}//{}{}/{}\n'.format(
            vpt2_info[3], vpt2_info[1], vpt2_info[2],
            vpt2_ref_info[3], vpt2_ref_info[1], vpt2_ref_info[2])
    if sp_str:
        chemkin_header_str += sp_str

    return chemkin_header_str


def run_ckin_poly(spc, spc_dct_i, pac99_poly_str):
    """ prepare chemkin header info and convert pac 99 format to chemkin format
    """

    hf_str = '! Hf(0 K) = {:.2f}, Hf(298 K) = {:.2f} kcal/mol\n'.format(
        float(spc_dct_i['Hfs'][0]), float(spc_dct_i['Hfs'][1]))
    ich = spc_dct_i['ich']
    formula = automol.inchi.formula(ich)
    atom_dict = thermo.util.get_atom_counts_dict(formula)
    chemkin_poly_str = thermo.nasapoly.convert_pac_to_chemkin(
        spc, atom_dict, hf_str, pac99_poly_str)
    print('\nCHEMKIN Polynomial:')
    print(chemkin_poly_str)
    return chemkin_poly_str


def write_nasa_file(spc_dct_i, ckin_path, nasa_path, chemkin_poly_str):
    """ write out the nasa polynomials
    """
    ich = spc_dct_i['ich']
    formula = automol.inchi.formula(ich)
    with open(os.path.join(nasa_path, formula+'.ckin'), 'w') as nasa_file:
        nasa_file.write(chemkin_poly_str)
    with open(os.path.join(ckin_path, formula+'.ckin'), 'w') as nasa_file:
        nasa_file.write(chemkin_poly_str)
