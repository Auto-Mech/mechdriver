"""
  Uses information to assess what model type is appropriate
"""

import automol
from mechlib.amech_io import printer as ioprinter


def nonrigid_rotations(spc_mod_dct_i):
    """ dtermine if a nonrigid rotation model is specified and further
        information is needed from the filesystem
    """
    rot_model = spc_mod_dct_i['rot']['mod']
    return bool(rot_model == 'vpt2')


def nonrigid_tors(spc_mod_dct_i, rotors):
    """ dtermine if a nonrigid torsional model is specified and further
        information is needed from the filesystem
    """
    vib_model = spc_mod_dct_i['vib']['mod']
    tors_model = spc_mod_dct_i['tors']['mod']
    has_tors = bool(any(rotors))
    tors_hr_model = bool(
        tors_model in ('1dhr', '1dhrf', '1dhrfa', 'mdhr', 'mdhrv'))
    tau_hr_model = bool(tors_model == 'tau' and vib_model != 'vib')
    # diatomic model?
    return has_tors and (tors_hr_model or tau_hr_model)


def anharm_vib(spc_mod_dct_i):
    """ a
    """
    vib_model = spc_mod_dct_i['vib']['mod']
    return bool(vib_model == 'vpt2')


def tau_pf(spc_mod_dct_i):
    """ determine if pf is done with tau
    """
    tors_model = spc_mod_dct_i['tors']['mod']
    return bool(tors_model == 'tau')


def scale_1d(spc_mod_dct_i):
    """ determine if we need to scale the potential
    """
    ioprinter.debug_message(
        'tors model in scale set', spc_mod_dct_i['tors']['mod'])
    return bool(spc_mod_dct_i['tors']['mod'] == '1dhrfa')


def scale_tors_pot(spc_mod_dct_i, to_scale):
    """ determine if we need to scale the potential
    """
    onedhr_model = bool('1dhr' in spc_mod_dct_i['tors']['mod'])
    return bool(onedhr_model and to_scale)


def vib_tau(spc_mod_dct_i):
    """ determine if vibrations are treated via tau sampling
    """
    vib_model = spc_mod_dct_i['vib']['mod']
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


def need_fake_wells(tsclass, well_model):
    """ Return boolean to see if fake wells are needed
    """
    if well_model == 'fake':
        abst_rxn = bool('abstraction' in tsclass)
        # addn_rxn = bool('addition' in tsclass)
        subs_rxn = bool('substitution' in tsclass)
        need = bool(abst_rxn or subs_rxn)
    else:
        need = False

    return need


def var_radrad(tsclass):
    """ Return boolean to see if fake wells are needed
    """
    rad_rad = 'radical-radical' in tsclass
    low_spin = 'high' not in tsclass

    return bool(rad_rad and low_spin)


def treat_tunnel(ts_mod, ts_class):
    """ decide to treat tunneling
    """

    treat = True

    ts_sadpt, ts_nobar = ts_mod['sadpt'], ts_mod['nobar']
    tunnel_model = ts_mod['tunnel']
    radrad = 'radical-radical' in ts_class
    if tunnel_model is not None:
        if radrad:
            if ts_nobar in ('pst', 'rpvtst', 'vrctst'):
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


def is_abstraction(spc_dct):
    """ Check if a PES consists of a single abstraction reaction.
    """
    _abstraction = False

    ts_names = tuple(name for name in spc_dct.keys() if 'ts_' in name)
    if len(ts_names) == 1:
        ts_dct = spc_dct[ts_names[0]]
        if 'abstraction' in ts_dct['zrxn'].class_:
            _abstraction = True

    return _abstraction
