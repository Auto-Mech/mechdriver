"""
  Write polrate file
"""

import os
from ioformat import build_mako_str
from polyrate_io._util import pt_format


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def potential_file(rct1_info, sadpt_info, mep_infos,
                   rct2_info=None, prd1_info=None, prd2_info=None):
    """ write the potential file
    """

    # Write the reactants
    rct1_hess, rct1_vval, rct1_sval = rct1_info
    rct_str = pt_format(
        'r140', rct1_hess, 'vvalue', rct1_vval,
        slabel='svalue', sval=rct1_sval, geo=None, grad=None)
    if rct2_info is not None:
        rct2_hess, rct2_vval, rct2_sval = rct2_info
        rct_str += pt_format(
            'r240', rct2_hess, 'vvalue', rct2_vval,
            slabel='svalue', sval=rct2_sval, geo=None, grad=None)

    # Write the product
    if prd1_info is not None:
        prd1_hess, prd1_vval, prd1_sval = prd1_info
        prd_str = pt_format(
            'p140', prd1_hess, 'vvalue', prd1_vval,
            slabel='svalue', sval=prd1_sval, geo=None, grad=None)
    if prd2_info is not None:
        prd2_hess, prd2_vval, prd2_sval = prd2_info
        prd_str += pt_format(
            'p240', prd2_hess, 'vvalue', prd2_vval,
            slabel='svalue', sval=prd2_sval, geo=None, grad=None)

    # Write the saddle point
    sadpt_hess, sadpt_vval = sadpt_info
    sadpt_str = pt_format(
        'saddle40', sadpt_hess, 'vvalue', sadpt_vval,
        slabel=None, sval=None, geo=None, grad=None)

    # Write the points along the MEP
    mep_str = ''
    for mep_info in mep_infos:
        hess, vval, sval, geo, grad = mep_info
        mep_str += pt_format(
            'point40', hess, 'vmep', vval,
            slabel='smep', sval=sval, geo=geo, grad=grad)
        mep_str += '\n'

    npts = len(mep_infos)

    dat_str = (
        rct_str + '\n' +
        # prd_str + '\n' +
        sadpt_str + '\n' +
        mep_str + '\n'
    )

    pot40_keys = {
        'title': 'Polyrate',
        'npts': npts,
        'dat_str': dat_str

    }

    return build_mako_str(
        template_file_name='pot40.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys=pot40_keys,
        remove_whitespace=False)
