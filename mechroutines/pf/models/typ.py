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
    tors_model = spc_mod_dct_i['vib']['mod']
    has_tors = bool(any(rotors))
    tors_hr_model = bool(
        tors_model in ('1dhr', '1dhrf', '1dhrfa', 'mdhr', 'mdhrv'))
    tau_hr_model = bool(tors_model == 'tau' and vib_model != 'vib')

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
        'tors model in scale set', spc_mod_dct_i['mod']['tors'])
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


def pst_ts(rxn_class, ts_sadpt, ts_nobarrier):
    """ Return boolean to see if fake wells are needed
    """

    pst = False
    if not automol.par.is_radrad(rxn_class):
        if ts_sadpt == 'pst':
            pst = True
    else:
        if ts_nobarrier == 'pst':
            pst = True

    return pst


def need_fake_wells(rxn_class, well_model):
    """ Return boolean to see if fake wells are needed
    """
    if well_model == 'fake':
        need = any(rtyp in automol.par.typ(rxn_class)
                   for rtyp in ('abstraction', 'substitution'))
    else:
        need = False

    return need


def treat_tunnel(ts_mod, rxn_class):
    """ decide to treat tunneling
    """

    treat = True

    ts_sadpt, ts_nobar = ts_mod['sadpt'], ts_mod['nobar']
    tunnel_model = ts_mod['tunnel']
    if tunnel_model is not None:
        if automol.par.is_radrad(rxn_class):
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
        rxn_class = spc_dct[ts_names[0]]['class']
        if 'abstraction' in automol.par.typ(rxn_class):
            _abstraction = True

    return _abstraction
