"""
  Uses information to assess what model type is appropriate
"""

import automol


def nonrigid_rotations(pf_models):
    """ dtermine if a nonrigid rotation model is specified and further
        information is needed from the filesystem
    """
    rot_model = pf_models['rot']
    return bool(rot_model == 'vpt2')


def nonrigid_tors(pf_models, rotors):
    """ dtermine if a nonrigid torsional model is specified and further
        information is needed from the filesystem
    """
    vib_model, tors_model = pf_models['vib'], pf_models['tors']
    has_tors = bool(any(rotors))
    tors_hr_model = bool(
        tors_model in ('1dhr', '1dhrf', '1dhrfa', 'mdhr', 'mdhrv'))
    tau_hr_model = bool(tors_model == 'tau' and vib_model != 'vib')
    # diatomic model?
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


def pst_ts(tsclass, ts_sadpt, ts_nobarrier):
    """ Return boolean to see if fake wells are needed
    """

    pst = False
    if not var_radrad(tsclass):
        if ts_sadpt == 'pst':
            pst = True
    else:
        if ts_nobarrier == 'pst':
            pst = True

    return pst


def need_fake_wells(tsclass):
    """ Return boolean to see if fake wells are needed
    """
    abst_rxn = bool('abstraction' in tsclass)
    # addn_rxn = bool('addition' in tsclass)
    subs_rxn = bool('substitution' in tsclass)
    return bool(abst_rxn or subs_rxn)
    # return bool(abst_rxn or addn_rxn or subs_rxn)


def var_radrad(tsclass):
    """ Return boolean to see if fake wells are needed
    """
    rad_rad = 'radical radical' in tsclass
    low_spin = 'high spin' not in tsclass
    addn_rxn = 'addition' in tsclass
    return bool(rad_rad and low_spin and addn_rxn)


def treat_tunnel(tunnel_model, ts_sadpt, ts_nobarrier, radrad):
    """ decide to treat tunneling
    """
    treat = True
    if tunnel_model != 'none':
        if radrad:
            if ts_nobarrier in ('pst', 'vrctst'):
                treat = False
        else:
            if ts_sadpt == ('pst', 'vrctst'):
                treat = False

    return treat


def is_atom(spc_dct_i):
    """ Check if species is an atom
    """
    geo = automol.inchi.geometry(spc_dct_i['inchi'])
    return automol.geom.is_atom(geo)
