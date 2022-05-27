"""
  Handle building the energy transfer section
"""

from phydat import phycon
import automol
import mess_io
from mechanalyzer.inf import spc as sinfo
from mechlib.amech_io import printer as ioprinter


# FUNCTIONS TO WRITE THE ENERGY TRANSFER  STRINGS
def make_energy_transfer_strs(tgt_info, bath_info, etrans_dct):
    """ Writes the energy down (`Edown`) and collision frequency (`ColldFreq`)
        section strings for the MESS input.

        Function and subfunctions examine the etrans_dct to decide how to
        set the energy transfer parameters written into the string:
        mass, lj params, and alpha params.

        If user-specified options in etrans_dct, those will be used. If
        'estimate' given in the etrans_dct. Then the function will use the
        bath and target spc_info objects to estimate the values using Jasper
        formulae.
    """

    # Determine all of the energy transfer parameters
    ioprinter.info_message(
        '- Determining the masses...', newline=1)
    mass1, mass2 = mass_params(
        tgt_info, bath_info, etrans_dct)

    ioprinter.info_message(
        '- Determining the Lennard-Jones model parameters...', newline=1)
    sig, eps = lj_params(
        tgt_info, bath_info, etrans_dct)
    sig1, eps1, sig2, eps2 = sig, eps, sig, eps
    # ^ Get two sets of A+B params. Are written in MESS and then uncombined

    ioprinter.info_message(
        '- Determining the energy-down transfer model parameters...',
        newline=1)
    exp_factor, exp_power, exp_cutoff = edown_params(
        tgt_info, bath_info, etrans_dct,
        ljpar=(sig1, eps1, mass1, mass2))

    # Write the Energy Transfer section string
    if all(val is not None
           for val in (sig1, sig2, eps1, eps2, exp_factor, exp_power)):
        edown_str = mess_io.writer.energy_down(
            exp_factor, exp_power, exp_cutoff)
        collid_freq_str = mess_io.writer.collision_frequency(
            eps1, eps2, sig1, sig2, mass1, mass2)
    else:
        # set defaults?
        edown_str, collid_freq_str = None, None

    return edown_str, collid_freq_str


# FUNCTIONS TO SET ALL OF THE PARAMETERS FOR THE MESS FILE
def mass_params(tgt_info, bath_info, etrans_dct):
    """ Determine the mass parameters used to define the energy transfer model
        for target-bath collisions in master equation simulations.

        Function first assesses if a user has supplied mass parameters
        by way of etrans_dct. If masses have not been provided, then the masses
        will be simply calculated using the geometries of the target and
        bath species.

        :param tgt_info: spc_info object for target species
        :type tgt_info: mechanalyzer.inf.spc object
        :param bath_info: spc_info object for bath species
        :type bath_info: mechanalyzer.inf.spc object
        :param etrans_dct: energy transfer options for mass, lj, alpha
        :type: etrans_dct: dict[str:tuple/str]
    """

    # First try and get the mass for spc or PES (global_mass)
    mass = etrans_dct.get('mass', None)
    if mass is None:
        mass = etrans_dct.get('global_mass', None)

    # Now set the values
    if mass is not None:
        ioprinter.info_message('  - Using the user input values...')
        [mass1, mass2] = mass
    else:
        ioprinter.info_message('  - Obtaining masses from geometries...')
        mass1 = sum(automol.geom.masses(automol.inchi.geometry(tgt_info[0])))
        mass2 = sum(automol.geom.masses(automol.inchi.geometry(bath_info[0])))

    return mass1, mass2


def lj_params(tgt_info, bath_info, etrans_dct):
    """ Determine the Lennard-Jones (LJ) epsilon and sigma parameters used to
        define the energy transfer model for target-bath collisions in
        master equation simulations.

        Function first assesses if a user has supplied LJ parameters
        by way of etrans_dct. If params have not been provided, then the params
        will be simply calculated estimating them using the Jasper formulae,
        which requires the InChI strings of the target and bath.

        :param tgt_info: spc_info object for target species
        :type tgt_info: mechanalyzer.inf.spc object
        :param bath_info: spc_info object for bath species
        :type bath_info: mechanalyzer.inf.spc object
        :param etrans_dct: energy transfer options for mass, lj, alpha
        :type: etrans_dct: dict[str:tuple/str]
    """

    sig, eps = None, None

    # First try and get the ljpar for spc or PES (global_ljpar)
    ljpar = etrans_dct.get('lj', 'estimate')
    if ljpar == 'estimate':
        ljpar = etrans_dct.get('global_lj', 'estimate')

    # Now set the values
    if ljpar is not None:

        if isinstance(ljpar, list):

            ioprinter.info_message('  - Using the user input values...')
            [sig, eps] = ljpar

        elif ljpar == 'estimate':

            ioprinter.info_message('- Estimating the parameters...')
            tgt_ich, bath_ich = tgt_info[0], bath_info[0]
            model = automol.etrans.estimate.determine_collision_model_series(
                tgt_ich, bath_ich, 'lj')
            n_heavy = automol.geom.atom_count(
                automol.inchi.geometry(tgt_ich), 'H', match=False)
            ioprinter.info_message(
                '    - Series to use for estimation for estimation: '
                f' {model}\n'
                f'    - Heavy atom count: {n_heavy}')

            sig, eps = automol.etrans.estimate.lennard_jones_params(
                n_heavy, model)
    else:
        sig, eps = None, None

    return sig, eps


def edown_params(tgt_info, bath_info, etrans_dct, ljpar=None):
    """ Determine the average collisional energy (alpha) parameter used to
        define the energy transfer model for target-bath collisions in
        master equation simulations.

        Function first assesses if a user has supplied alpha parameter
        by way of etrans_dct. If params have not been provided, then the param
        will be simply calculated estimating them using the Jasper formulae,
        which requires the InChI strings of the target and bath, as well as
        masses and LJ params.

        :param tgt_info: spc_info object for target species
        :type tgt_info: mechanalyzer.inf.spc object
        :param bath_info: spc_info object for bath species
        :type bath_info: mechanalyzer.inf.spc object
        :param etrans_dct: energy transfer options for mass, lj, alpha
        :type: etrans_dct: dict[str:tuple/str]
    """

    efactor, epower, ecutoff = None, None, None

    # First try and get the edown for spc or PES (global_edown)
    edown = etrans_dct.get('edown', 'estimate')
    if edown == 'estimate':
        edown = etrans_dct.get('global_edown', 'estimate')

    # Now get the values
    if edown is not None:

        if isinstance(edown, list):

            ioprinter.info_message('  - Using the user input values...')
            [efactor, epower, ecutoff] = edown

        elif edown == 'estimate':

            assert ljpar is not None

            ioprinter.info_message('  - Estimating the parameters...')

            tgt_geo = automol.inchi.geometry(tgt_info[0])
            model = automol.etrans.estimate.determine_collision_model_series(
                tgt_info[0], bath_info[0], 'alpha')
            if model is not None:
                sig, eps, mass1, mass2 = ljpar
                n_eff = automol.etrans.estimate.effective_rotor_count(tgt_geo)
                ioprinter.info_message(
                    '    - Series to use for estimation for estimation:'
                    f' {model}\n'
                    f'    - Found effective rotor count: {n_eff:.2f}\n'
                    '    - Using following LJ parameters for '
                    'collisional frequency and alpha calculation:\n'
                    f'       eps={eps*phycon.EH2WAVEN:.2f} cm-1, '
                    f'sigma={sig*phycon.BOHR2ANG:.2f} Ang,\n'
                    f'       mass1={mass1:.2f} amu, mass2={mass2:.2f} amu')
                efactor, epower = automol.etrans.estimate.alpha(
                    n_eff, eps, sig, mass1, mass2, model,
                    empirical_factor=2.0)
                ecutoff = 15.0
            else:
                efactor, epower, ecutoff = None, None, None
                ioprinter.warning_message(
                    f'Cannot calculate Zalpha for {model}')

    return efactor, epower, ecutoff


# FUNCTION TO HANDLE BUILDING ETRANS OBJECTS
def etrans_dct_for_species(spc_dct_i, pes_mod_dct_i):
    """  Build an energy transfer dictionary for species from a spc dct
    """
    return {
        'bath': pes_mod_dct_i['energy_transfer']['bath'],
        'mass': spc_dct_i.get('mass', None),
        'lj': spc_dct_i.get('lj', 'estimate'),
        'edown': spc_dct_i.get('edown', 'estimate')
    }


def set_etrans_well(rxn_lst, spc_dct):
    """ Determine the species that can represent energy transfer
        for the entire PES.
    """

    _, (reacs, prods) = rxn_lst[0]
    if len(reacs) == 1:
        well_dct = spc_dct[reacs[0]]
        well_name = reacs[0]
    elif len(prods) == 1:
        well_dct = spc_dct[prods[0]]
        well_name = prods[0]
    else:
        rct1_dct = spc_dct[reacs[0]]
        rct2_dct = spc_dct[reacs[1]]
        rct1_count = automol.geom.count(
            automol.inchi.geometry(rct1_dct['inchi']))
        rct2_count = automol.geom.count(
            automol.inchi.geometry(rct2_dct['inchi']))
        if rct1_count > rct2_count:
            well_dct = rct1_dct
            well_name = reacs[0]
        else:
            well_dct = rct2_dct
            well_name = reacs[1]

    well_info = sinfo.from_dct(well_dct)
    ioprinter.info_message(
        f'  - Using {well_name} for global collision target for PES')

    return well_info


def set_bath(spc_dct, etrans_dct):
    """ Build nfo object for the bath
    """

    # Try to obtain bath set by the user, otherwise use N2
    bath_name = etrans_dct.get('bath', None)
    bath_dct = spc_dct.get(bath_name, None)
    if bath_dct is not None:
        bath_info = sinfo.from_dct(bath_dct)
        ioprinter.info_message(f'  - Using bath {bath_name} input by user')
    else:
        bath_info = ['InChI=1S/Ar', 0, 1]
        ioprinter.info_message(
            '  - No bath provided, using Argon bath as default')

    return bath_info
