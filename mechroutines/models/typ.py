"""
  Uses information to assess what model type is appropriate
"""

import automol
from mechlib.amech_io import printer as ioprinter


def nonrigid_rotations(spc_mod_dct_i):
    """ Determine if the rotational partition function for a certain
        species should be calculated according to some non-rigid model.

        This determination solely relies on whether has specified the
        use of a non-rigid model for the species.

        :param spc_mod_dct_i: species partition function models
        :type spc_mod_dct_i: dict[str: str]
        :rtype: bool
    """
    rot_model = spc_mod_dct_i['rot']['mod']
    return bool(rot_model == 'vpt2')


def nonrigid_tors(spc_mod_dct_i, rotors):
    """ Determine if internal rotation portions of the vibrational
        partition function both can and should be treated with some
        non-rigid torsion model.

        This is is determined by (1) if the user has requested the use
        of a non-rigid model, and (2) whether internal rotors exist for
        the given species.

        :param spc_mod_dct_i: species partition function models
        :type spc_mod_dct_i: dict[str: str]
        :rtype: bool
    """
    vib_model = spc_mod_dct_i['vib']['mod']
    tors_model = spc_mod_dct_i['tors']['mod']
    has_tors = bool(rotors is not None)
    tors_hr_model = bool(
        tors_model in ('1dhr', '1dhrf', '1dhrfa', 'mdhr', 'mdhrv',
                       'tau-1dhr', 'tau-1dhrf', 'tau-1dhrfa')
    )
    tau_hr_model = bool('tau' in tors_model and vib_model != 'vib')
    if has_tors:
        if any([len(pot)<1 for pot in automol.rotor.potentials(rotors, flat=True)]):
            print('WARNING: empty potential will crash MESS so using rigid model instead')
            has_tors = False
        
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
    return bool(tors_model in ('tau', 'tau-1dhr', 'tau-1dhrf', 'tau-1dhrfa'))


def scale_1d(spc_mod_dct_i):
    """ determine if we need to scale the potential
    """
    ioprinter.debug_message(
        'tors model in scale set', spc_mod_dct_i['tors']['mod'])
    return (
        bool(spc_mod_dct_i['tors']['mod'] in ('1dhrfa', '1dhrf', '1dhr'))
        and spc_mod_dct_i['tors']['scale'] == 'on')


def scale_tors_pot(spc_mod_dct_i, to_scale):
    """ determine if we need to scale the potential
    """
    onedhr_model = bool('1dhr' in spc_mod_dct_i['tors']['mod'])
    return bool(onedhr_model and to_scale)


def squash_tors_pot(spc_mod_dct_i):
    """ determine if we need to scale the potential
    """
    ioprinter.debug_message(
        'tors model in scale set', spc_mod_dct_i['tors']['mod'])
    return bool(spc_mod_dct_i['tors']['mod'] in ('1dhrfa', 'tau-1dhrfa'))


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
    """ Determine if master equation treatments of a reaction channel
        necessitate the generation and inclusion of fake van der Waals wells
        at the entrance- and exit-channel.

        This is is determined by (1) if the user has requested the use
        of fake-wells, and (2) whether the reaction is bimol???

        :param rxn_class:
        :type rxn_class:
        :param well_model:
        :type well_model:
        :rtype: bool
    """
    return (well_model == 'fake' and automol.par.need_wells(rxn_class))


def treat_tunnel(ts_mod, rxn_class, ts_inf_dct=None):
    """ Determine if master equation treatments of a reaction channel
        should include treatments of quantum tunneling treatments.

        This is is determined by (1) if the user has requested some tunneling
        treatment, and the transition state corresponds to a saddle point.

        :param ts_mod:
        :type ts_mod:
        :rtype: bool
    """

    treat = True

    ts_sadpt, ts_nobar = ts_mod['sadpt'], ts_mod['nobar']
    tunnel_model = ts_mod['tunnel']
    if ts_inf_dct:
        if ts_inf_dct['writer'] == 'pst_block':
            treat = False
    if tunnel_model is not None:
        if automol.par.is_radrad(rxn_class):
            if ts_nobar in ('pst', 'rpvtst', 'vrctst'):
                treat = False
        else:
            if ts_sadpt == ('pst', 'vrctst'):
                treat = False
    return treat


def is_atom(spc_dct_i):
    """ Determine if the species is an atom by using the information
        provided in its data dictionary of the overall species dictionary.

        :param spc_dct_i:
        :type spc_dct_i:
        :rtype: bool
    """
    # return automol.geom.is_atom(automol.chi.geometry(spc_dct_i['inchi']))
    return sum(automol.chi.formula(spc_dct_i['inchi']).values()) < 2


def is_abstraction_pes(spc_dct, rxn_lst, pes_idx):
    """ Assess whether a potential energy surface consists of a single
        abstraction reaction.

        will break???
    """

    _abstraction = False

    if len(rxn_lst) == 1:
        chnl_idx, _ = rxn_lst[0]
        tsname = f'ts_{pes_idx+1:g}_{chnl_idx+1:g}_0'

        rxn_class = spc_dct[tsname]['class']
        if (automol.par.typ(rxn_class) ==
           automol.par.ReactionClass.Typ.HYDROGEN_ABSTRACTION):
            _abstraction = True

    return _abstraction
