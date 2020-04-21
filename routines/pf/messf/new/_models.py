"""
  Build the species string based on the model
"""

import elstruct
import automol
import autofile
from routines.es.scan import hr_prep
from routines.pf.messf import _rot as rot
from routines.pf.messf import _tors as tors
from routines.pf.messf import _sym as sym
from routines.pf.messf import _vib as vib
from routines.pf.messf import _tau as tau
# from routines.pf.messf import _vpt2 as vpt2
from lib.phydat import phycon
from lib.filesystem import orb as fsorb
from lib.filesystem import minc as fsmin
from lib.filesystem import build as fbuild
from lib.filesystem import inf as finf


def read_filesys_for_spc(spc_dct_i, rxn, spc_model, pf_levels,
                         save_prefix, saddle=False):
    """ Pull all of the neccessary information from the filesystem for a species
    """

    # Unpack the models and levels
    rot_model, tors_model, vib_model, sym_model = spc_model

    # Set filesys
    cnf_save_fs, cnf_save_path, cnf_save_locs, thy_save_path = _cnf_filesys(
        spc_dct_i, rxn, pf_levels, save_prefix, saddle=saddle, level='harm')

    # Initialize all of the elemetnts of the inf dct
    geom, sym_factor, freqs, imag, elec_levels = None, None, None, None, None
    has_tors, hr_str, mdhr_mess_str, mdhr_dat_str = False, '', '', ''
    xmat, rovib_coups, rot_dists = None, None, None

    # PULL ALL INFORMATION FROM THE FILESYS AND WRITE THE MESS SPECIES BLOCK

    # Rotation: Geometry
    geom = rot.read_geom(cnf_save_fs, cnf_save_locs)

    # Rotation: Anharmonicity
    if nonrigid_rotations(rot_model):
        vpt2_save_fs, _, vpt2_save_locs, _ = _cnf_filesys(
            spc_dct_i, rxn, pf_levels, save_prefix,
            saddle=saddle, level='vpt2')
        rovib_coups, rot_dists = rot.read_rotational_values(
            vpt2_save_fs, vpt2_save_locs)

    # Torsion Info: Needed for Proper symmetry and vibration determination
    if nonrigid_tors(vib_model, tors_model, tors_names):
        mess_hr_str, proj_hr_str, mdhr_dat_str = tors.proc_mess_strings()
        # mdhr_dat_file_name = '{}_mdhr.dat'.format(spc[0])

    # Symmetry Factor
    sym_factor = sym.symmetry_factor(
        sym_model, spc_dct_i, spc_info, dist_names,
        saddle, frm_bnd_key, brk_bnd_key, tors_names,
        tors_cnf_save_fs, tors_cnf_save_locs,
        sym_cnf_save_fs, sym_cnf_save_locs)
    if nonrigid_tors(vib_model, tors_model, tors_names):
        sym_factor = sym.tors_reduced_sym_factor(
            cnf_save_locs, cnf_save_fs, saddle=False)

    # Vibrations: Frequencies
    if nonrigid_tors(vib_model, tors_model, tors_names):
        freqs, imag_freq, zpe = tors.tors_freqs_zpve()  # Could maybe split up
    else:
        freqs, imag = vib.read_harmonic_freqs(
            geom, cnf_save_fs, cnf_save_locs, saddle=saddle)
    
    # Vibrations: Anharmonicity
    if anharm_vib:
        xmat = vib.anharmonicity()

    # Elec levels
    elec_levels = spc_dct_i['elec_levels']

    # Create info dictionary
    keys = ['geom', 'sym_factor', 'freqs', 'imag', 'elec_levels',
            'has_tors', 'hr_str', 'mdhr_mess_str', 'mdhr_dat_str',
            'xmat', 'rovib_coups', 'rot_dists']
    vals = [geom, sym_factor, freqs, imag, elec_levels,
            has_tors, hr_str, mdhr_mess_str, mdhr_dat_str,
            xmat, rovib_coups, rot_dists]
    inf_dct = dict(zip(keys, vals))

    return inf_dct


# VRCTST
def read_filesys_for_flux(ts_dct, rxn, pf_levels, save_prefix):
    """ Grab the flux file from the filesystem
    """

    # Set filesys
    ts_save_fs, ts_save_path = _ts_filesys(
        spc_dct, rxn, pf_levels, save_prefix, level='harm')

    # Read the flux file string
    locs = []
    flux_str = ts_save_fs[-1].file.flux.read(locs)

    return flux_str


# VTST
def read_filesys_for_rpvtst(rpath_vals, sadpt=True):
    """ Pull all of the neccessary information from the filesystem for a species
    """
    # Set filesys
    if sadpt:
        _, cnf_save_path, _, _ = _cnf_filesys(
            spc_dct_i, rxn, pf_levels, save_prefix,
            saddle=saddle, level='harm')
        scn_save_fs, scn_locs, save_paths = _scn_filesys(
            cnf_save_path, run_tors_names)
    else:
        ts_save_fs, ts_save_path = _ts_filesys(
            spc_dct, rxn, pf_levels, save_prefix, level='harm')
        scn_save_fs, scn_locs, save_paths = _scn_filesys(
            ts_save_path, run_tors_names)

    # Loop over scan filesystem and pull out the values
    inf_dct_lst = []
    for locs in scn_locs:

        # Check if to pull info
        if locs not in rpath_vals:
            continue

        # Get geometry, energy, vibrational freqs, and zpe
        if scn_save_fs[-1].file.geometry.exists(locs):
            geom = scn_save_fs[-1].file.geometry.read(locs)
        else:
            print('no geom')
            continue
        if scn_save_fs[-1].file.energy.exists(locs):
            sp_save_fs = autofile.fs.single_point(scn_save_path)
            sp_level = fsorb.mod_orb_restrict(ts_info, ene_thy_level)
            if sp_save_fs[-1].file.energy.exists(sp_level[1:4]):
                ene = sp_save_fs[-1].file.energy.read(sp_level[1:4])
            else:
                print('no energy')
                continue
        else:
            print('no energy')
            continue
        if scn_save_fs[-1].file.hessian.exists(locs):
            proj_rotors_str = ''
            hess = scn_save_fs[-1].file.hessian.read(locs)
            scn_save_path = scn_save_fs[-1].path(locs)
            freqs, _, _ = vib.projrot_freqs_1(
                geom, hess,
                proj_rotors_str,
                scn_save_path, pot=False, saddle=True)
            zpe = sum(freqs)*phycon.WAVEN2KCAL/2.
        else:
            print('no hessian')
            continue

        # Get the relative energy
        zero_ene = ''

        # Create info dictionary and append to lst
        sym_factor = 1.0
        elec_levels = ts_dct['elec_levels']
        keys = ['geom', 'sym_factor', 'freqs', 'elec_levels', 'zero_ene']
        vals = [geom, sym_factor, freqs, elec_levels, zero_ene]
        inf_dct_lst.append(dict(zip(keys, vals)))

    return inf_dct_lst


# TAU
def read_filesys_for_tau(spc_dct_i, spc_model, pf_levels,
                         save_prefix, rxn=(), saddle=False):
    """ Read the filesystem to get information for TAU
    """

    # Use model to determine whether to read grads and hessians
    _, vib_model, _ = spc_model
    if vib_model != 'tau':
        read_gradient, read_hessian = False, False
        freqs = ()
    else:
        read_gradient, read_hessian = True, True
        freqs = ()

    # Set the filesystem
    cnf_save_fs, _, cnf_save_locs, _ = _cnf_filesys(
        spc_dct_i, rxn, pf_levels, save_prefix, saddle=saddle, level='harm')
    tau_save_fs = autofile.fs.tau(save_prefix)

    # Read the reference geometry and energy
    if cnf_save_locs:
        ref_geom = cnf_save_fs[-1].file.geometry.read(cnf_save_locs)
        min_ene = cnf_save_fs[-1].file.energy.read(cnf_save_locs)

    # Set the ground and reference energy to set values for now
    ground_ene = -0.02
    reference_ene = 0.00

    # Write the flux mode str
    flux_mode_str = tors.write_flux_str()

    # Read the geom, ene, grad, and hessian for each sample
    samp_geoms, samp_enes, samp_grads, samp_hessians = [], [], [], []
    for locs in tau_save_fs[-1].existing():

        geo = tau_save_fs[-1].file.geometry.read(locs)
        geo_str = autofile.file.write.geometry(geo)
        samp_geoms.append(geo_str)

        ene = tau_save_fs[-1].file.energy.read(locs)
        rel_ene = (ene - min_ene) * phycon.EH2KCAL
        ene_str = autofile.file.write.energy(rel_ene)
        samp_enes.append(ene_str)

        if read_gradient:
            grad = tau_save_fs[-1].file.gradient.read(locs)
            grad_str = autofile.file.write.gradient(grad)
            samp_grads.append(grad_str)

        if read_hessian:
            hess = tau_save_fs[-1].file.hessian.read(locs)
            hess_str = autofile.file.write.hessian(hess)
            samp_hessians.append(hess_str)

    # Create info dictionary
    keys = ['ref_geom', 'elec_levels', 'freqs', 'flux_mode_str',
            'ground_ene', 'reference_ene',
            'samp_geoms', 'samp_enes', 'samp_grads', 'samp_hessians']
    vals = [ref_geom, spc_dct_i['elec_levels'], freqs, flux_mode_str,
            ground_ene, reference_ene,
            samp_geoms, samp_enes, samp_grads, samp_hessians]
    inf_dct = dict(zip(keys, vals))

    return inf_dct


# Filesystem object setters
def _cnf_filesys(spc_dct, rxn, pf_levels, save_prefix,
                 saddle=False, level='harm'):
    """ Set needed conformer filesys objects
    """

    if level == 'harm':
        thy_info = pf_levels[2]
    elif level == 'vpt2':
        thy_info = pf_levels[2]

    # Set the filesystem objects
    spc_info = (spc_dct['ich'], spc_dct['chg'], spc_dct['mul'])
    mod_thy_info = fsorb.mod_orb_restrict(spc_info, thy_info)
    if not saddle:
        _, thy_save_path = fbuild.spc_thy_fs_from_root(
            save_prefix, spc_info, mod_thy_info)
    else:
        rxn_info = finf.rxn_info(
            rxn['reacs'], rxn['prods'], spc_dct)
        _, thy_save_path = fbuild.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_thy_info)
        _, thy_save_path = fbuild.ts_fs_from_thy(thy_save_path)

    # Get the cnf filesys needed everything based off geo+freq (also vpt2)
    cnf_save_fs, cnf_save_locs = fbuild.cnf_fs_from_prefix(
        thy_save_path, cnf='min')
    cnf_save_paths = fbuild.cnf_paths_from_locs(
        cnf_save_fs, cnf_save_locs)

    return cnf_save_fs, cnf_save_paths[0], cnf_save_locs, thy_save_path


def _scn_filesys(cnf_save_path, run_tors_names):
    """ Set needed conformer filesys objects
    """

    # Get a list of the other tors coords to freeze and set the filesystem
    frz_all_tors = es_keyword_dct['frz_all_tors']
    if frz_all_tors:
        scn_save_fs = autofile.fs.cscan(cnf_save_path)
        scn_locs = fbuild.cscn_locs_from_fs(scn_save_fs, run_tors_names)
    else:
        scn_save_fs = autofile.fs.scan(cnf_save_path)
        scn_locs = fbuild.scn_locs_from_fs(scn_save_fs, run_tors_names)

    return scn_save_fs, scn_locs


def _ts_filesys(spc_dct, rxn, pf_levels, save_prefix, level='harm'):
    """ Set needed conformer filesys objects
    """

    if level == 'harm':
        thy_info = pf_levels[2]
    elif level == 'vpt2':
        thy_info = pf_levels[2]

    # Set the filesystem objects
    spc_info = (spc_dct['ich'], spc_dct['chg'], spc_dct['mul'])
    mod_thy_info = fsorb.mod_orb_restrict(spc_info, thy_info)
    rxn_info = finf.rxn_info(
        rxn['reacs'], rxn['prods'], spc_dct)
    _, thy_save_path = fbuild.rxn_thy_fs_from_root(
        save_prefix, rxn_info, mod_thy_info)
    ts_save_fs, ts_save_path = fbuild.ts_fs_from_thy(thy_save_path)

    return ts_save_fs, ts_save_path


# Series of checks to determine what information is needed to be obtained
def nonrigid_rotations(rot_model):
    """ dtermine if a nonrigid rotation model is specified and further
        information is needed from the filesystem
    """
    return bool(rot_model in ('vpt2'))


def nonrigid_tors(vib_model, tors_model, tors_names):
    """ dtermine if a nonrigid torsional model is specified and further
        information is needed from the filesystem
    """
    has_tors = bool(tors_names)
    tors_hr_model = bool(tors_model in ('1dhr', '1dhrf', 'mdhr', 'mdhrv'))
    tau_hr_model = bool(tors_model == 'tau' and vib_model != 'vib')
    # diatomic model
    return has_tors and (tors_hr_model or tau_hr_model)


def anharm_vib(vib_model):
    """ a
    """
    return bool(vib_model == 'vpt2')


def tau_pf(tors_model):
    """ determine if pf is done with tau
    """
    return bool(tors_model == 'tau')


def vib_tau(vib_model):
    """ determine if vibrations are treated via tau sampling
    """
    return bool(vib_model == 'tau')
