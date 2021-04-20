"""
Random functions that are needed in drivers and routines
"""

from mechroutines.pf.models import typ
from mechlib.amech_io import parser
from mechlib.amech_io import printer as ioprinter
from mechanalyzer.inf import rxn as rinfo


def set_pf_info(model_dct, thy_dct, chn_model, ref_model):
    """ set pf information
    """
    chn_pf_levels = parser.model.pf_level_info(
        model_dct[chn_model]['es'], thy_dct)
    chn_pf_models = parser.model.pf_model_info(
        model_dct[chn_model]['pf'])
    ref_pf_levels = parser.model.pf_level_info(
        model_dct[ref_model]['es'], thy_dct)
    ref_pf_models = parser.model.pf_model_info(
        model_dct[ref_model]['pf'])

    chn_pf_models['ref_scheme'] = (model_dct[chn_model]['options']['ref_scheme']
        if 'ref_scheme' in model_dct[ref_model]['options'] else 'none')
    ref_pf_models['ref_scheme'] = (model_dct[ref_model]['options']['ref_scheme']
        if 'ref_scheme' in model_dct[ref_model]['options'] else 'none')
    chn_pf_models['ref_enes'] = (model_dct[chn_model]['options']['ref_enes']
        if 'ref_enes' in model_dct[ref_model]['options'] else 'none')
    ref_pf_models['ref_enes'] = (model_dct[ref_model]['options']['ref_enes']
        if 'ref_enes' in model_dct[ref_model]['options'] else 'none')

    return chn_pf_levels, chn_pf_models, ref_pf_levels, ref_pf_models


def set_ts_cls_info(spc_dct, model_dct, tsname, chn_model):
    """ figure out various information needed to establish reaction class
    """
    sub_tsname = tsname + '_0'
    ts_class = spc_dct[sub_tsname]['zrxn'].class_
    ts_sadpt = model_dct[chn_model]['pf']['ts_sadpt']
    ts_nobarrier = model_dct[chn_model]['pf']['ts_barrierless']
    tunnel_model = model_dct[chn_model]['pf']['tunnel']
    _ = typ.var_radrad(ts_class)

    radrad = rinfo.radrad(spc_dct[sub_tsname]['rxn_info'])
    if radrad:
        ts_class = 'radical radical ' + ts_class

    return ts_class, ts_sadpt, ts_nobarrier, tunnel_model, radrad


def make_rxn_str(rlst, prepend=''):
    """ convert list to string
    """
    return prepend + '+'.join(rlst)


# Printing functions
def print_pf_info(pf_models, pf_levels, chn_model, ref_ene_lvl):
    """ printing stuff
    """

    ioprinter.info_message(
        'Model name for Channel: {}'.format(chn_model), newline=1)

    ioprinter.info_message('Partition Functions Treatments:')
    for key, val in pf_models.items():
        model_str = '  {} = {}'.format(key, val)
        ioprinter.info_message(model_str)

    ioprinter.info_message(
        'Electronic Structure Levels for all Required Data:')
    level_str = ''
    for key, val in pf_levels.items():
        if key == 'ene':
            if val[0] != ref_ene_lvl:
                level_str = '  {} = {}'.format(key, val[0])
                level_str += '  * differs from reference; will calc shifts'
        if key == 'tors':
            level_str = '  tors_sp = {}\n'.format(val[0][0])
            level_str += ' tors_scn = {}'.format(val[0][1])
        else:
            level_str = '  {} = {}'.format(key, val[0])
        ioprinter.info_message(level_str)
