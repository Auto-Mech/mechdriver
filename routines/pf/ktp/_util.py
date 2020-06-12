"""
  Utilitiy functions
"""

from lib.amech_io import parser


def set_pf_info(model_dct, thy_dct, chn_model, ref_model):
    """ set pf information
    """
    chn_pf_levels = parser.model.set_es_model_info(
        model_dct[chn_model]['es'], thy_dct)
    chn_pf_models = parser.model.set_pf_model_info(
        model_dct[chn_model]['pf'])
    ref_pf_levels = parser.model.set_es_model_info(
        model_dct[ref_model]['es'], thy_dct)
    ref_pf_models = parser.model.set_pf_model_info(
        model_dct[ref_model]['pf'])

    return chn_pf_levels, chn_pf_models, ref_pf_levels, ref_pf_models


def set_ts_cls_info(spc_dct, model_dct, tsname, chn_model):
    """ figure out various information needed to establish reaction class
    """
    ts_class = spc_dct[tsname]['class']
    ts_sadpt = model_dct[chn_model]['pf']['ts_sadpt']
    ts_nobarrier = model_dct[chn_model]['pf']['ts_barrierless']
    tunnel_model = model_dct[chn_model]['pf']['tunnel']

    return ts_class, ts_sadpt, ts_nobarrier, tunnel_model


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


def make_rxn_str(rlst, prepend=''):
    """ convert list to string
    """
    return prepend + '+'.join(rlst)


# def treat_tunnel(tun_model, ts_sadpt, ts_nobarrier, _var_radrad(ts_class)):
def treat_tunnel():
    """ decide to treat tunneling
    """
    return True
