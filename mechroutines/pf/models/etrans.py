"""
  Handle building the energy transfer section
"""

import automol
import mess_io
from mechanalyzer.inf import spc as sinfo
from mechroutines.pf.models import _eff as eff
from mechlib.amech_io import printer as ioprinter


# FUNCTIONS TO WRITE THE ENERGY TRANSFER  STRINGS
def make_energy_transfer_strs(well_info, bath_info, etrans_dct):
    """ Makes the standard header and energy transfer sections
        for MESS input file
    """

    # Determine all of the energy transfer parameters
    ioprinter.info_message(
        '- Determining the masses...', newline=1)
    mass1, mass2 = mass_params(
        well_info, bath_info, etrans_dct)

    print(etrans_dct)

    ioprinter.info_message(
        '- Determining the Lennard-Jones model parameters...', newline=1)
    sig1, eps1, sig2, eps2 = lj_params(
        well_info, bath_info, etrans_dct)

    print(sig1, eps1, sig2, eps2)
    ioprinter.info_message(
        '- Determining the energy-down transfer model parameters...', 
        newline=1)
    exp_factor, exp_power, exp_cutoff = edown_params(
        well_info, bath_info, etrans_dct,
        ljpar=(sig1, eps1, mass1, mass2))
    print(exp_factor, exp_power, exp_cutoff)

    # Write the Energy Transfer section string
    if all(val is not None
           for val in (sig1, sig2, eps1, eps2, exp_factor, exp_power)):
        edown_str = mess_io.writer.energy_down(
            exp_factor, exp_power, exp_cutoff)
        collid_freq_str = mess_io.writer.collision_frequency(
            eps1, eps2, sig1, sig2, mass1, mass2)
    else:
        edown_str, collid_freq_str = None, None

    return edown_str, collid_freq_str


# FUNCTIONS TO SET ALL OF THE PARAMETERS FOR THE MESS FILE
def mass_params(well_info, bath_info, etrans_dct):
    """ Build the mass parameters
    """

    mass = etrans_dct.get('mass', None)
    if mass is not None:

        ioprinter.info_message('  - Using the user input values...')
        [mass1, mass2] = mass

    else:

        ioprinter.info_message('  - Obtaining masses from geometries...')
        ioprinter.info_message('well', well_info[0])
        geo = automol.inchi.geometry(well_info[0])
        mass1 = sum(automol.geom.masses(geo))

        geo = automol.inchi.geometry(bath_info[0])
        mass2 = sum(automol.geom.masses(geo))

    return mass1, mass2


def lj_params(well_info, bath_info, etrans_dct):
    """ Build the lennard-jones parameters
    """

    sig1, eps1, sig2, eps2 = None, None, None, None

    ljpar = etrans_dct.get('lj', None)

    if ljpar is not None:

        if isinstance(ljpar, list):

            ioprinter.info_message('  - Using the user input values...')
            [sig1, eps1, sig2, eps2] = ljpar

        elif ljpar == 'estimate':

            ioprinter.info_message('- Estimating the parameters...')
            well_ich = well_info[0]
            well_geo = automol.inchi.geometry(well_ich)
            params = eff.estimate_viable(well_ich, well_geo, bath_info)
            if params is not None:
                bath_model, tgt_model = params
                ioprinter.info_message(
                    '    - Series to use for estimation for estimation...')
                ioprinter.info_message(
                    '      Bath: {}, Target: {} '.format(
                        bath_model, tgt_model))

                ioprinter.info_message(
                    '    - Effective atom numbers for estimation...')
                n_heavy = automol.geom.atom_count(well_geo, 'H', match=False)
                ioprinter.info_message(
                    '      N_heavy: ', n_heavy)
                sig, eps = eff.lj_sig_eps(n_heavy, bath_model, tgt_model)
                sig1, eps1, sig2, eps2 = sig, eps, sig, eps

        elif ljpar == 'read':

            ioprinter.info_message('- Reading the filesystem...')
            ljs_lvl = etrans_dct.get('ljlvl', None)
            if ljs_lvl is not None:
                # Get the levels into theory objects
                pf_filesystems = 0
                sig1, eps1, sig2, eps2 = _read_lj(pf_filesystems)

    else:
        sig1, eps1, sig2, eps2 = None, None, None, None

    return sig1, eps1, sig2, eps2


def edown_params(well_info, bath_info, etrans_dct, ljpar=None):
    """ Energy down model parameters
    """

    efactor, epower, ecutoff = None, None, None

    edown = etrans_dct.get('edown', None)

    if edown is not None:

        if isinstance(edown, list):

            ioprinter.info_message('  - Using the user input values...')
            [efactor, epower, ecutoff] = edown

        elif edown == 'estimate':

            assert ljpar is not None

            ioprinter.info_message('  - Estimating the parameters...')
            well_ich = well_info[0]
            well_geo = automol.inchi.geometry(well_ich)
            params = eff.estimate_viable(well_ich, well_geo, bath_info)
            if params is not None:
                bath_model, tgt_model = params
                ioprinter.info_message(
                    '    - Series to use for estimation for estimation...')
                ioprinter.info_message(
                    '      Bath: {}, Target: {} '.format(bath_model, tgt_model))

                ioprinter.info_message(
                    '    - Calculating the LJ collisional frequencies...')
                sig, eps, mass1, mass2 = ljpar
                zlj_dct = eff.lj_collision_frequency(
                    sig, eps, mass1, mass2,
                    temps=(300., 1000., 2000.))
                ioprinter.info_message(
                    '    - Effective atom numbers for estimation...')
                n_eff = eff.calc_n_eff(well_geo)
                ioprinter.info_message(
                    '      N_eff: ', n_eff)
                efactor, epower = eff.alpha(n_eff, zlj_dct,
                                            bath_model, tgt_model)
                ecutoff = 15.0

        elif edown == 'read':

            ioprinter.info_message('  - Reading the filesystem...')
            edownlvl = etrans_dct.get('edownlvl', None)
            if edownlvl is not None:
                # NEED: Get the levels into theory objects
                pf_filesystems = 0
                efactor = _read_alpha(pf_filesystems)

    return efactor, epower, ecutoff


# FILESYS READING
def _read_alpha(pf_filesystems):
    """ filesys
    """

    etrans_fs, etrans_locs = 0, 0
    # etrans_fs, etrans_locs = filesys.models.set_etrans_fs()

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

    etrans_fs, etrans_locs = 0, 0
    # etrans_save_fs, etrans_locs = filesys.models.set_etrans_fs()

    # Read the epsilon and sigma values
    if etrans_fs[-1].file.alpha.exists(etrans_locs):
        alpha = etrans_fs[-1].file.alpha.read(etrans_locs)
    else:
        alpha = None

    return alpha


# FUNCTION TO HANDLE BUILDING ETRANS OBJECTS
def build_etrans_dct(spc_dct_i):
    """  Build an energy transfer dict from a spc dct
    """

    etrans_dct = {}

    mass = spc_dct_i.get('mass', None)
    if mass is not None:
        etrans_dct['mass'] = mass

    ljpar = spc_dct_i.get('lj', None)
    if ljpar is not None:
        etrans_dct['lj'] = ljpar

    edown = spc_dct_i.get('edown', None)
    if edown is not None:
        etrans_dct['edown'] = edown

    return etrans_dct


def set_etrans_well(rxn_lst, spc_dct):
    """ Build info object for reference well on PES
    """

    well_dct = None

    reacs = rxn_lst[0]['reacs']
    prods = rxn_lst[0]['prods']
    if len(reacs) == 1:
        well_dct = spc_dct[reacs[0]]
    elif len(prods) == 1:
        well_dct = spc_dct[prods[0]]
    else:
        rct1_dct = spc_dct[reacs[0]]
        rct2_dct = spc_dct[reacs[1]]
        rct1_count = automol.geom.count(
            automol.inchi.geometry(rct1_dct['inchi']))
        rct2_count = automol.geom.count(
            automol.inchi.geometry(rct2_dct['inchi']))
        if rct1_count > rct2_count:
            well_dct = rct1_dct
        else:
            well_dct = rct2_dct

    well_info = sinfo.from_dct(well_dct)

    return well_info


def set_bath(spc_dct, etrans_dct):
    """ Build nfo object for the bath
    """

    # Try to obtain bath set by the user, otherwise use N2
    bath_name = etrans_dct.get('bath', None)
    bath_dct = spc_dct.get(bath_name, None)
    if bath_dct is not None:
        bath_info = sinfo.from_dct(bath_dct)
        ioprinter.info_message(
            '  - Using bath {} input by user'.format(bath_name))
    else:
        bath_info = ['InChI=1S/Ar', 0, 1]
        ioprinter.info_message(
            '  - No bath provided, using Argon as bath')

    return bath_info
