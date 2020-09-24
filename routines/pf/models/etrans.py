"""
  Handle building the energy transfer section
"""

import automol
from routines.pf.models import _eff as eff
from lib.filesys import inf as finf


# FUNCTION TO HANDLE THE BATH SETTING
def set_well(rxn_lst, spc_dct):
    """ Build info object for reference well on PES
    """
    reac = rxn_lst[0]['reacs'][0]
    reac_info = finf.get_spc_info(spc_dct[reac])

    return reac_info


def set_bath(spc_dct, etrans_dct):
    """ Build nfo object for the bath
    """

    # Try to obtain bath set by the user, otherwise use N2
    bath_name = etrans_dct.get('bath', None)
    bath_dct = spc_dct.get(bath_name, None)
    if bath_dct is not None:
        bath_info = finf.get_spc_info(bath_dct)
        print('  - Using bath {} input by user'.format(bath_name))
    else:
        bath_info = ['InChI=1S/N2/c1-2', 0, 1]
        print('  - No bath provided, using N2 as bath')

    return bath_info


# FUNCTIONS TO SET ALL OF THE PARAMETERS FOR THE MESS FILE
def edown_params(well_info, bath_info, spc_dct, etrans_dct):
    """ Energy down model parameters
    """
    
    edown = etrans_dct.get('edown', None)
    if isinstance(edown, list):

        print('  - Using the user input values...')
        [efactor, epower, ecutoff] = edown
    
    elif edown == 'estimate':

        print('  - Estimating the parameters...')
        well_geo = automol.inchi.geometry(well_info[0])
        params = eff.estimate_viable(well_geo, bath_info)
        if params is not None:
            bath_model, tgt_model = params
            print('    - Series to use for estimation for estimation...')
            print('      Bath: {}, Target: {} '.format(bath_model, tgt_model))

            print('    - Determining the effective atom numbers for estimation...')
            n_eff = eff.calc_n_eff(well_geo)
            n_heavy = automol.geom.atom_count(well_geo, 'H', match=False)
            print('      N_eff: ', n_eff)
            print('      N_heavy: ', n_heavy)
            efactor, epower = eff.alpha(n_eff, n_heavy, bath_model, tgt_model)
            ecutoff = 15.0

    elif edown == 'read':

        print('  - Reading the filesystem...')
        edownlvl = etrans_dct.get('edownlvl', None)
        if edownlvl is not None:
            # NEED: Get the levels into theory objects
            efactor = _read_alpha(pf_filesystems)
        
    return efactor, epower, ecutoff


def lj_params(well_info, bath_info, spc_dct, etrans_dct):
    """ Build the lennard-jones parameters
    """

    lj = etrans_dct.get('lj', None)
    if isinstance(lj, list):

        print('  - Using the user input values...')
        [sig1, eps1, sig2, eps2] = lj
    
    elif lj == 'estimate':

        print('- Estimating the parameters...')
        well_geo = automol.inchi.geometry(well_info[0])
        params = eff.estimate_viable(well_geo, bath_info)
        if params is not None:
            bath_model, tgt_model = params
            print('    - Series to use for estimation for estimation...')
            print('      Bath: {}, Target: {} '.format(bath_model, tgt_model))
            
            print('    - Determining the effective atom numbers for estimation...')
            n_heavy = automol.geom.atom_count(well_geo, 'H', match=False)
            print('      N_heavy: ', n_heavy)
            sig, eps = eff.lj(n_heavy, bath_model, tgt_model)
            sig1, eps1, sig2, eps2 = sig, eps, sig, eps
    
    elif lj == 'read':

        print('- Reading the filesystem...')
        ljs_lvl = etrans_dct.get('ljlvl', None)
        if ljs_lvl is not None:
            # Get the levels into theory objects
            sig1, eps1, sig2, eps2 = _read_lj(pf_filesystems)
        
    return sig1, eps1, sig2, eps2


def mass_params(well_info, bath_info, spc_dct, etrans_dct):
    """ Build the mass parameters
    """

    mass = etrans_dct.get('mass', None)
    if mass is not None:

        print('  - Using the user input values...')
        [mass1, mass2] = mass

    else:

        print('  - Obtaining masses from geometries...')
        geo = automol.inchi.geometry(well_info[0])
        mass1 = sum(automol.geom.masses(geo))

        geo = automol.inchi.geometry(bath_info[0])
        mass2 = sum(automol.geom.masses(geo))

    return mass1, mass2


# FILESYS READING
def _read_alpha(pf_filesystems):
    """ filesys
    """

    etrans_save_fs, etrans_locs = filesys.set_etrans_fs()

    # Read the epsilon and sigma values
    if etrans_fs[-1].file.epsilon.exists(etrans_locs):
        eps = etrans_fs[-1].file.epsilon.read(etrans_locs)
    else:
        eps = None
    if etrans_fs[-1].file.sigma.exists(etrans_locs):
        sig = etrans_fs[-1].file.sigma.read(etrans_locs)
    else:
        sig = None 

    return eps, sig


def _read_lj(pf_filesystems):
    """ filesys
    """

    etrans_save_fs, etrans_locs = filesys.set_etrans_fs()

    # Read the epsilon and sigma values
    if etrans_fs[-1].file.alpha.exists(etrans_locs):
        alpha = etrans_fs[-1].file.alpha.read(etrans_locs)
    else:
        alpha = None

    return alpha
