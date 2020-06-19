"""
  Read the save filesystem for all of the required information specified by
    (1) the models specified for partition function and
    (2) the electronic structure levels
  in order to write portions of MESS strings for species and reaction paths
  and calculate electronic and zero-point vibrational energies.
"""

import autofile
from routines.pf.messf import ene
from routines.pf.messf import _rot as rot
from routines.pf.messf import _tors as tors
from routines.pf.messf import _sym as sym
from routines.pf.messf import _vib as vib
from routines.pf.messf import _util as util
from lib.filesys import inf as finf
from lib.filesys import mincnf


def atm_data(spc_dct_i,
             chn_pf_models, chn_pf_levels, ref_pf_models, ref_pf_levels,
             run_prefix, save_prefix):
    """ Pull all neccessary info for the atom
    """

    # Set up all the filesystem objects using models and levels
    pf_filesystems = build_pf_filesystems(
        spc_dct_i, chn_pf_levels, run_prefix, save_prefix, False)

    print('\nObtaining the geometry...')
    geom = rot.read_geom(pf_filesystems)
    
    print('\nObtaining the electronic energy...')
    ene_chnlvl = ene.read_energy(
        spc_dct_i, pf_filesystems, chn_pf_models, chn_pf_levels,
        read_ene=True, read_zpe=False)

    ene_reflvl = None
    _, _ = ref_pf_models, ref_pf_levels

    # Create info dictionary
    inf_dct = {
        'geom': geom,
        'sym_factor': 1.0,
        'freqs': [],
        'mess_hr_str': '',
        'mass': util.atom_mass(spc_dct_i),
        'elec_levels': spc_dct_i['elec_levels'],
        'ene_chnlvl': ene_chnlvl,
        'ene_reflvl': ene_reflvl
    }

    return inf_dct


def mol_data(spc_dct_i,
             chn_pf_models, chn_pf_levels, ref_pf_models, ref_pf_levels,
             run_prefix, save_prefix, saddle=False, tors_wgeo=False):
    """ Pull all of the neccessary information from the filesystem for a species
    """

    # Initialize all of the elements of the inf dct
    geom, sym_factor, freqs, imag, elec_levels = None, None, None, None, None
    mess_hr_str, mdhr_dats = '', []
    xmat, rovib_coups, rot_dists = None, None, None

    # Set up all the filesystem objects using models and levels
    pf_filesystems = build_pf_filesystems(
        spc_dct_i, chn_pf_levels, run_prefix, save_prefix, saddle)

    # Set information for transition states
    frm_bnd_keys, brk_bnd_keys = util.get_bnd_keys(pf_filesystems, saddle)
    rxn_class = util.set_rxn_class(spc_dct_i, saddle)

    # Obtain rotor information used to determine new information
    print('\nPreparing internal rotor info building partition functions...')
    rotor_names, rotor_grids, rotor_syms, const_dct, ref_ene = tors.rotor_info(
        spc_dct_i, pf_filesystems, chn_pf_models,
        frm_bnd_key=frm_bnd_keys, brk_bnd_key=brk_bnd_keys)

    if nonrigid_tors(chn_pf_models, rotor_names):
        mess_hr_str, prot_hr_str, mdhr_dats, rotor_syms = tors.make_hr_strings(
            rotor_names, rotor_grids, rotor_syms, const_dct,
            ref_ene, pf_filesystems, chn_pf_models,
            rxn_class, frm_bnd_keys,
            saddle=saddle, tors_wgeo=tors_wgeo)

    # Obtain rotation partition function information
    print('\nObtaining info for rotation partition function...')
    geom = rot.read_geom(pf_filesystems)

    if nonrigid_rotations(chn_pf_models):
        rovib_coups, rot_dists = rot.read_rotational_values(pf_filesystems)

    # Obtain vibration partition function information
    print('\nObtaining the vibrational frequencies and zpves...')
    if nonrigid_tors(chn_pf_models, rotor_names):
        freqs, imag, zpe = vib.tors_projected_freqs_zpe(
            pf_filesystems, mess_hr_str, prot_hr_str, saddle=saddle)
    else:
        freqs, imag, zpe = vib.read_harmonic_freqs(
            pf_filesystems, saddle=saddle)

    if anharm_vib(chn_pf_models):
        xmat = vib.read_anharmon_matrix(pf_filesystems)

    # Obtain symmetry factor
    print('\nDetermining the symmetry factor...')
    sym_factor = sym.symmetry_factor()
    #     sym_model, spc_dct_i, spc_info, dist_names,
    #     saddle, frm_bnd_key, brk_bnd_key, rotor_names,
    #     cnf_save_fs, cnf_save_locs, saddle)

    if nonrigid_tors(chn_pf_models, rotor_names):
        sym_factor = sym.tors_reduced_sym_factor(
            sym_factor, rotor_syms)

    # Obtain electronic energy levels
    elec_levels = spc_dct_i['elec_levels']

    # Obtain energy levels
    print('\nObtaining the electronic energy...')
    chn_ene = ene.read_energy(
        spc_dct_i, pf_filesystems, chn_pf_models, chn_pf_levels,
        read_ene=True, read_zpe=False)
    # print('mod ene', chn_ene, zpe)
    ene_chnlvl = chn_ene + zpe

    ene_reflvl = None
    _, _ = ref_pf_models, ref_pf_levels
    # if chn_model == ref_model:
    #     ene_reflvl = ene_chnlvl
    # else:
    #     ene_reflvl = get_fs_ene_zpe(spc_dct, prod,
    #                                 thy_dct, model_dct, model,
    #                                 save_prefix, saddle=False,
    #                                 read_ene=True, read_zpe=True)

    # Create info dictionary
    keys = ['geom', 'sym_factor', 'freqs', 'imag', 'elec_levels',
            'mess_hr_str', 'mdhr_dats',
            'xmat', 'rovib_coups', 'rot_dists',
            'ene_chnlvl', 'ene_reflvl']
    vals = [geom, sym_factor, freqs, imag, elec_levels,
            mess_hr_str, mdhr_dats,
            xmat, rovib_coups, rot_dists,
            ene_chnlvl, ene_reflvl]
    inf_dct = dict(zip(keys, vals))

    return inf_dct


# VRCTST
# def flux_data(ts_dct, rxn, pf_levels, save_prefix):
#     """ Grab the flux file from the filesystem
#     """
#
#     # Set filesys
#     ts_save_fs, ts_save_path = _ts_filesys(
#         spc_dct, rxn, pf_levels, save_prefix, level='harm')
#
#     # Read the flux file string
#     locs = []
#     flux_str = ts_save_fs[-1].file.flux.read(locs)
#
#     return flux_str


# VTST
# def rpvtst_data(rpath_vals, sadpt=True):
#     """ Pull all of the neccessary information from the
#         filesystem for a species
#     """
#     # Set filesys
#     if sadpt:
#         _, cnf_save_path, _, _ = _cnf_filesys(
#             spc_dct_i, rxn, pf_levels, save_prefix,
#             saddle=saddle, level='harm')
#         scn_save_fs, scn_locs, save_paths = _scn_filesys(
#             cnf_save_path, run_tors_names)
#     else:
#         ts_save_fs, ts_save_path = _ts_filesys(
#             spc_dct, rxn, pf_levels, save_prefix, level='harm')
#         scn_save_fs, scn_locs, save_paths = _scn_filesys(
#             ts_save_path, run_tors_names)
#
#     # Loop over scan filesystem and pull out the values
#     inf_dct_lst = []
#     for locs in scn_locs:
#
#         # Check if to pull info
#         if locs not in rpath_vals:
#             continue
#
#         # Get geometry, energy, vibrational freqs, and zpe
#         if scn_save_fs[-1].file.geometry.exists(locs):
#             geom = scn_save_fs[-1].file.geometry.read(locs)
#         else:
#             print('no geom')
#             continue
#         if scn_save_fs[-1].file.energy.exists(locs):
#             sp_save_fs = autofile.fs.single_point(scn_save_path)
#             sp_level = fsorb.mod_orb_restrict(ts_info, ene_thy_level)
#             if sp_save_fs[-1].file.energy.exists(sp_level[1:4]):
#                 ene = sp_save_fs[-1].file.energy.read(sp_level[1:4])
#             else:
#                 print('no energy')
#                 continue
#         else:
#             print('no energy')
#             continue
#         if scn_save_fs[-1].file.hessian.exists(locs):
#             proj_rotors_str = ''
#             hess = scn_save_fs[-1].file.hessian.read(locs)
#             scn_save_path = scn_save_fs[-1].path(locs)
#             freqs, _, _ = vib.projrot_freqs_1(
#                 geom, hess,
#                 proj_rotors_str,
#                 scn_save_path, pot=False, saddle=True)
#             zpe = sum(freqs)*phycon.WAVEN2KCAL/2.
#         else:
#             print('no hessian')
#             continue
#
#         # Get the relative energy
#         zero_ene = ''
#
#         # Create info dictionary and append to lst
#         sym_factor = 1.0
#         elec_levels = ts_dct['elec_levels']
#         keys = ['geom', 'sym_factor', 'freqs', 'elec_levels', 'zero_ene']
#         vals = [geom, sym_factor, freqs, elec_levels, zero_ene]
#         inf_dct_lst.append(dict(zip(keys, vals)))
#
#     return inf_dct_lst


# TAU
# def tau_data(spc_dct_i, spc_name,
#              chn_pf_models, chn_pf_levels, ref_pf_models, ref_pf_levels,
#              run_prefix, save_prefix, saddle=False):
#     """ Read the filesystem to get information for TAU
#     """
#
#     # Use model to determine whether to read grads and hessians
#     _, vib_model, _ = pf_models
#     if vib_model != 'tau':
#         read_gradient, read_hessian = False, False
#         freqs = ()
#     else:
#         read_gradient, read_hessian = True, True
#         freqs = ()
#
#     # Set the filesystem
#     cnf_save_fs, _, cnf_save_locs, _ = _cnf_filesys(
#         spc_dct_i, rxn, pf_levels, save_prefix, saddle=saddle, level='harm')
#     tau_save_fs = autofile.fs.tau(save_prefix)
#
#     # Read the reference geometry and energy
#     if cnf_save_locs:
#         ref_geom = cnf_save_fs[-1].file.geometry.read(cnf_save_locs)
#         min_ene = cnf_save_fs[-1].file.energy.read(cnf_save_locs)
#
#     # Set the ground and reference energy to set values for now
#     ground_ene = -0.02
#     reference_ene = 0.00
#
#     # Write the flux mode str
#     flux_mode_str = tors.write_flux_str()
#
#     # Read the geom, ene, grad, and hessian for each sample
#     samp_geoms, samp_enes, samp_grads, samp_hessians = [], [], [], []
#     for locs in tau_save_fs[-1].existing():
#
#         geo = tau_save_fs[-1].file.geometry.read(locs)
#         geo_str = autofile.file.write.geometry(geo)
#         samp_geoms.append(geo_str)
#
#         ene = tau_save_fs[-1].file.energy.read(locs)
#         rel_ene = (ene - min_ene) * phycon.EH2KCAL
#         ene_str = autofile.file.write.energy(rel_ene)
#         samp_enes.append(ene_str)
#
#         if read_gradient:
#             grad = tau_save_fs[-1].file.gradient.read(locs)
#             grad_str = autofile.file.write.gradient(grad)
#             samp_grads.append(grad_str)
#
#         if read_hessian:
#             hess = tau_save_fs[-1].file.hessian.read(locs)
#             hess_str = autofile.file.write.hessian(hess)
#             samp_hessians.append(hess_str)
#
#     # Create info dictionary
#     keys = ['ref_geom', 'elec_levels', 'freqs', 'flux_mode_str',
#             'ground_ene', 'reference_ene',
#             'samp_geoms', 'samp_enes', 'samp_grads', 'samp_hessians']
#     vals = [ref_geom, spc_dct_i['elec_levels'], freqs, flux_mode_str,
#             ground_ene, reference_ene,
#             samp_geoms, samp_enes, samp_grads, samp_hessians]
#     inf_dct = dict(zip(keys, vals))
#
#     return inf_dct


# Filesystem object creators
def build_pf_filesystems(spc_dct_i, pf_levels,
                         run_prefix, save_prefix, saddle):
    """ Create various filesystems needed
    """

    pf_filesystems = {}
    # print('pf_levels', pf_levels)
    pf_filesystems['harm'] = set_model_filesys(
        spc_dct_i, pf_levels['harm'][1], run_prefix, save_prefix, saddle)
    if pf_levels['sym']:
        pf_filesystems['sym'] = set_model_filesys(
            spc_dct_i, pf_levels['sym'][1], run_prefix, save_prefix, saddle)
    if pf_levels['tors']:
        pf_filesystems['tors'] = set_model_filesys(
            spc_dct_i, pf_levels['tors'][1][0],
            run_prefix, save_prefix, saddle)
    if pf_levels['vpt2']:
        pf_filesystems['vpt2'] = set_model_filesys(
            spc_dct_i, pf_levels['vpt2'][1], run_prefix, save_prefix, saddle)

    return pf_filesystems


def set_model_filesys(spc_dct_i, level, run_prefix, save_prefix, saddle):
    """ Gets filesystem objects for torsional calculations
    """

    # Set the spc_info
    spc_info = finf.get_spc_info(spc_dct_i)

    # Set some path stuff
    if saddle:
        save_path = spc_dct_i['rxn_fs'][3]
        run_path = spc_dct_i['rxn_fs'][2]
    else:
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_fs[-1].create(spc_info)
        save_path = spc_save_fs[-1].path(spc_info)
        spc_run_fs = autofile.fs.species(run_prefix)
        spc_run_fs[-1].create(spc_info)
        run_path = spc_run_fs[-1].path(spc_info)

    # Set theory filesystem used throughout
    thy_save_fs = autofile.fs.theory(save_path)
    thy_run_fs = autofile.fs.theory(run_path)

    # Set the level for the model
    levelp = finf.modify_orb_restrict(spc_info, level)

    # Get the save fileystem path
    save_path = thy_save_fs[-1].path(levelp[1:4])
    run_path = thy_run_fs[-1].path(levelp[1:4])
    if saddle:
        save_fs = autofile.fs.transition_state(save_path)
        save_fs[0].create()
        save_path = save_fs[0].path()
        run_fs = autofile.fs.transition_state(run_path)
        run_fs[0].create()
        run_path = run_fs[0].path()

    # Get the fs object and the locs
    # print('save path', save_path)
    cnf_save_fs = autofile.fs.conformer(save_path)
    min_cnf_locs = mincnf.min_energy_conformer_locators(cnf_save_fs)

    # Get the save path for the conformers
    if min_cnf_locs:
        cnf_save_path = cnf_save_fs[-1].path(min_cnf_locs)
    else:
        cnf_save_path = ''

    return [cnf_save_fs, cnf_save_path, min_cnf_locs, save_path, run_path]


# Series of checks to determine what information is needed to be obtained
def nonrigid_rotations(pf_models):
    """ dtermine if a nonrigid rotation model is specified and further
        information is needed from the filesystem
    """
    rot_model = pf_models['rot']
    return bool(rot_model in ('vpt2'))


def nonrigid_tors(pf_models, tors_names):
    """ dtermine if a nonrigid torsional model is specified and further
        information is needed from the filesystem
    """
    vib_model, tors_model = pf_models['vib'], pf_models['tors']
    has_tors = bool(any(tors_names))
    tors_hr_model = bool(tors_model in ('1dhr', '1dhrf', 'mdhr', 'mdhrv'))
    tau_hr_model = bool(tors_model == 'tau' and vib_model != 'vib')
    # diatomic model
    return has_tors and (tors_hr_model or tau_hr_model)


def anharm_vib(pf_models):
    """ a
    """
    vib_model = pf_models['vib']
    return bool(vib_model == 'vpt2')


def tau_pf(pf_models):
    """ determine if pf is done with tau
    """
    tors_model = pf_models['tors']
    return bool(tors_model == 'tau')


def vib_tau(pf_models):
    """ determine if vibrations are treated via tau sampling
    """
    vib_model = pf_models['vib']
    return bool(vib_model == 'tau')
