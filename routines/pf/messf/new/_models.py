"""
  Build the species string based on the model
"""

import elstruct
import automol
from routines.es.scan import hr_prep
from routines.pf.messf import _tors as tors
from routines.pf.messf import _vib as vib
from routines.pf.messf import _tau as tau
# from routines.pf.messf import _vpt2 as vpt2
from lib.phydat import phycon


def read_filesys_for_spc(spc_dct, spc_info, spc_model,
                         pf_levels, save_prefix, saddle=False):
    """ Pull all of the neccessary information from the filesystem for a species
    """

    # Unpack the models and levels
    tors_model, vib_model, sym_model = spc_model

    # Set filesys
    cnf_save_fs, cnf_save_path, cnf_save_locs, thy_save_path = _cnf_filesys(
        spc_dct, pf_levels, level='harm')

    # Build the species string and get the imaginary frequency
    # Initialize various auxiliary MESS data strings and names
    mdhr_dat_str = ''
    mdhr_dat_file_name = '{}_mdhr.dat'.format(spc[0])
    # Initialize various variables and strings
    symf = sym_factor
    imag = 0.0
    freqs = []
    hr_str = ""
    xmat = ()

    # Default strings
    hind_ot_str, proj_rot_str = '', ''

    # Pull all information from the filesys and write the MESS species block
    if messfutil.is_atom(harm_min_cnf_locs, harm_cnf_save_fs):
        mass = messfutil.atom_mass(harm_min_cnf_locs, harm_cnf_save_fs)
        spc_str = mess_io.writer.atom(
            mass, elec_levels)
    else:
        # Geo (Rotational)
        if nonrigid_rotations(rot_model):
            rovib_coup, rot_dists = rot.read_rotational_values(
                vpt2_save_fs, vpt2_min_cnf_locs)
        else:
            rovib_coup, rot_dists = (), ()

        # Torsions
        if nonrigid_tors(vib_model, tors_model, tors_names):
            mess_hr_str, proj_hr_str, mdhr_dat_str = tors.proc_mess_strings()
        else:
            hind_rot_str, proj_rotors_str, mdhr_dat_str = '', '', ''

        # Symmetry
        sym_factor = sym.symmetry_factor(
            sym_model, spc_dct_i, spc_info, dist_names,
            saddle, frm_bnd_key, brk_bnd_key, tors_names,
            tors_cnf_save_fs, tors_min_cnf_locs,
            sym_cnf_save_fs, sym_min_cnf_locs)
        if nonrigid_tors(vib_model, tors_model, tors_names):
            sym_nums = tors.get_tors_sym_nums(
                spc_dct_i, tors_min_cnf_locs, tors_cnf_save_fs,
                frm_bnd_key, brk_bnd_key, saddle=False)
            for num in sym_nums:
                sym_factor /= num

        # Vibrations
        harm_geo, freqs, imag_freq = harm_freqs()
        freqs, imag_freq, zpe = tors.tors_freqs_zpve()  # Could maybe split up
        if anharm_vib:
            xmat = vib.anharmonicity()
        else:
            xmat = ()

    return spc_str, dat_str_dct, imag


# VTST
def read_filesys_for_rpvtst():
    """ Pull all of the neccessary information from the filesystem for a species
    """
    for idx in irc_idxs:

        # Set the filesystem locators for each grid point
        locs = [[dist_name], [idx]]
        print(scn_save_fs[-1].path(locs))

        # Get geometry, energy, vibrational freqs, and zpe
        if scn_save_fs[-1].file.geometry.exists(locs):
            geom = scn_save_fs[-1].file.geometry.read(locs)
        else:
            print('no geom')
            continue
        if scn_save_fs[-1].file.energy.exists(locs):
            if ene_thy_level == geo_thy_level:
                ene = scn_save_fs[-1].file.energy.read(locs)
                print('ene', ene)
            else:
                scn_save_path = scn_save_fs[-1].path(locs)
                sp_save_fs = autofile.fs.single_point(scn_save_path)
                sp_level = fsorb.mod_orb_restrict(ts_info, ene_thy_level)
                if sp_save_fs[-1].file.energy.exists(sp_level[1:4]):
                    ene = sp_save_fs[-1].file.energy.read(sp_level[1:4])
                    print('ene-high', ene)
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

# TAU
def read_filesys_for_tau():
    """ Read the filesystem to get information for TAU
    """

    # Read the filesystem for data string
    cnf_save_fs = autofile.fs.conformer(save_prefix)
    min_cnf_locs = fmin.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        ene_ref = cnf_save_fs[-1].file.energy.read(min_cnf_locs)
    tau_save_fs = autofile.fs.tau(save_prefix)

    # Read the geom, ene, grad, and hessian for each sample
    samp_geoms, samp_enes, samp_grads, samp_hessians = [], [], [], []
    for locs in tau_save_fs[-1].existing():

        geo = tau_save_fs[-1].file.geometry.read(locs)
        geo_str = autofile.file.write.geometry(geo)
        samp_geoms.append(geo_str)

        ene = tau_save_fs[-1].file.energy.read(locs)
        ene = (ene - ene_ref) * phycon.EH2KCAL
        ene_str = autofile.file.write.energy(ene)
        samp_enes.append(ene_str)

        if read_gradient:
            grad = tau_save_fs[-1].file.gradient.read(locs)
            grad_str = autofile.file.write.gradient(grad)
            samp_grads.append(grad_str)

        if read_hessian:
            hess = tau_save_fs[-1].file.hessian.read(locs)
            hess_str = autofile.file.write.hessian(hess)
            samp_hessians.append(hess_str)

    # Add the values to the dictionary
    inf_dct['samp_geoms'] = samp_geoms
    inf_dct['samp_enes'] = samp_enes         
    inf_dct['samp_grads'] = samp_grads       
    inf_dct['samp_hessians'] = samp_hessians          
              
    return inf_dct


# Filesystem object setters
def _cnf_filesys(spc_dct, pf_levels, level='harm'):
    """ """
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
            spc['reacs'], spc['prods'], spc_dct)
        _, thy_save_path = fbuild.rxn_thy_fs_from_root(
            save_prefix, rxn_info, mod_thy_info)
        _, thy_save_path = fbuild.ts_fs_from_thy(thy_save_path)

    # Get the cnf filesys needed everything based off geo+freq (also vpt2)
    cnf_save_fs, cnf_save_locs = fbuild.cnf_fs_from_prefix(
        thy_save_path, cnf='min')
    cnf_save_paths = fbuild.cnf_paths_from_locs(
        cnf_save_fs, cnf_save_locs)

    return cnf_save_fs, cnf_save_paths[0], cnf_save_locs, thy_save_path


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
