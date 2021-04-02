""" Library of parser functions for the model file
"""

import sys
import copy
import autoparse.find as apf
from ioformat import ptt
from mechanalyzer.inf import thy as tinfo
from mechlib.amech_io.parser.keywords import MODEL_PF_SUPPORTED_DCT
from mechlib.amech_io.parser.keywords import MODEL_PF_DEFAULT_DCT


MODEL_INP = 'inp/models.dat'


# FUNCTION TO READ IN A STRING FOR A SPECIFIC LEVEL #
def build_pes_model_keyword_dct(model_str):
    """ Build a dictionary for all the models keywords
    """

    # New grabbers
    pressures_str = _paren_block('pressures', model_str)
    rate_temps_str = _paren_block('rate_temps', model_str)
    thermo_temps_str = _paren_block('thermo_temps', model_str)
    # glob_etrans_str = _paren_block('etransfer', model_str)
    # rate_fit_str = _paren_block('rate_fit', model_str)
    # thermo_fit_str = _paren_block('thermo_fit', model_str)

    # Parse data
    pressures = ptt.keyword_lst(pressures_str)
    # OR pressures = ptt.build_vals_lst(pressures_str)
    rate_temps = ptt.keyword_lst(rate_temps_str)
    thermo_temps = ptt.keyword_lst(thermo_temps_str)
    # glob_etrans = ptt.keyword_dct(glob_etrans_str)
    # rate_fit = ptt.keyword_dct(rate_fit_str)
    # therm_fit = ptt.keyword_dct(therm_fit_str)

    return model_dct


def build_spc_model_keyword_dct(model_str):
    """ Build a dictionary for all the models keywords
    """
    # Grab the various sections required for each model
    pf_str = apf.first_capture(ptt.paren_section('pf'), model_str)
    es_str = apf.first_capture(ptt.paren_section('es'), model_str)

    # Get the dictionary for each section and check them
    pf_dct = ptt.build_keyword_dct(pf_str)
    es_dct = ptt.build_keyword_dct(es_str)

    # Combine dcts into single model dct
    model_dct = {}
    model_dct['pf'] = pf_dct
    model_dct['es'] = es_dct

    # Check the dct
    check_spc_model_dct(model_dct)

    return model_dct


def check_spc_model_dct(model_dct):
    """ Make sure the models dictionary keywords are all correct
    """

    # Check the pf dct
    pf_dct = model_dct['pf']
    for key, val in pf_dct.items():
        if key in MODEL_PF_SUPPORTED_DCT:
            if val not in MODEL_PF_SUPPORTED_DCT[key]:
                print('*ERROR: Value {}'.format(val),
                      'for Keyword {} not supported'.format(key))
                print('\nSupported keys for {}:'.format(key))
                for supp_val in MODEL_PF_SUPPORTED_DCT[key]:
                    print(supp_val)
                sys.exit()
        else:
            print('*ERROR: Keyword {} not supported'.format(key))
            print('\nSupported keys:')
            for supp_key in MODEL_PF_SUPPORTED_DCT:
                print(supp_key)
            sys.exit()

    # See if any pf model combinations were specified that are not supported
    check_model_combinations(pf_dct)

    # Check the es dct
    # es_dct = model_dct['es']

    # Check the vrctst dct
    # vrctst_dct = model_dct['vrctst']


def check_model_combinations(pf_dct):
    """ Check if a model combination is not implemented for PF routines
    """
    if pf_dct['vib'] == 'vpt2' and pf_dct['tors'] == '1dhr':
        print('*ERROR: VPT2 and 1DHR combination is not yet implemented')
        sys.exit()
    elif pf_dct['vib'] == 'vpt2' and pf_dct['tors'] == 'tau':
        print('*ERROR: VPT2 and TAU combination is not yet implemented')
        sys.exit()


def pf_model_info(pf_model):
    """ Set the PF model list based on the input
        Combine with es

        {'vib': {'model': '', 'geolvl_info': (), ...}}
    """
    rot_model = pf_model['rot'] if 'rot' in pf_model else 'rigid'
    tors_model = pf_model['tors'] if 'tors' in pf_model else 'rigid'
    vib_model = pf_model['vib'] if 'vib' in pf_model else 'harm'
    sym_model = pf_model['sym'] if 'sym' in pf_model else 'none'
    vpt2_model = pf_model['vpt2'] if 'vpt2' in pf_model else 'none'
    etrans_model = pf_model['etrans'] if 'etrans' in pf_model else 'none'

    # Set well models
    if 'wells' in pf_model:
        rwells_model = pf_model['wells']
        pwells_model = pf_model['wells']
    else:
        if 'rwells' in pf_model:
            rwells_model = pf_model['rwells']
        else:
            rwells_model = 'fake'
        if 'pwells' in pf_model:
            pwells_model = pf_model['pwells']
        else:
            pwells_model = 'fake'

    pf_models = {
        'rot': rot_model,
        'tors': tors_model,
        'vib': vib_model,
        'sym': sym_model,
        'vpt2': vpt2_model,
        'etrans': etrans_model,
        'rwells': rwells_model,
        'pwells': pwells_model,
    }

    return pf_models


def pf_level_info(es_model, thy_dct):
    """ Set the model info
    """
    # Read the ES models from the model dictionary
    geo_lvl = es_model['geo'] if 'geo' in es_model else None
    ene_lvl = es_model['ene'] if 'ene' in es_model else None
    harm_lvl = es_model['harm'] if 'harm' in es_model else None
    vpt2_lvl = es_model['vpt2'] if 'vpt2' in es_model else None
    sym_lvl = es_model['sym'] if 'sym' in es_model else None
    etrans_lvl = es_model['etrans'] if 'etrans' in es_model else None

    # Torsions and rxn paths which needs a reference for itself
    tors_lvl_sp = es_model['tors'][0] if 'tors' in es_model else None
    tors_lvl_scn = es_model['tors'][1] if 'tors' in es_model else None
    rpath_lvl_sp = es_model['rpath'][0] if 'rpath' in es_model else None
    rpath_lvl_scn = es_model['rpath'][1] if 'rpath' in es_model else None
    rpath_lvl_sp2 = es_model['rpath'][2] if 'rpath' in es_model else None

    # Set the theory info objects
    geo_thy_info = tinfo.from_dct(thy_dct.get(geo_lvl))
    harm_thy_info = tinfo.from_dct(thy_dct.get(harm_lvl))
    vpt2_thy_info = (tinfo.from_dct(thy_dct.get(vpt2_lvl))
                     if vpt2_lvl else None)
    sym_thy_info = (tinfo.from_dct(thy_dct.get(sym_lvl))
                    if sym_lvl else None)
    etrans_thy_info = (tinfo.from_dct(thy_dct.get(etrans_lvl))
                    if etrans_lvl else None)
    tors_sp_thy_info = (tinfo.from_dct(thy_dct.get(tors_lvl_sp))
                        if tors_lvl_sp else None)
    tors_scn_thy_info = (tinfo.from_dct(thy_dct.get(tors_lvl_scn))
                         if tors_lvl_scn else None)
    rpath_sp_thy_info = (tinfo.from_dct(thy_dct.get(rpath_lvl_sp))
                         if rpath_lvl_sp else None)
    rpath_scn_thy_info = (tinfo.from_dct(thy_dct.get(rpath_lvl_scn))
                          if rpath_lvl_scn else None)
    rpath_sp2_thy_info = (tinfo.from_dct(thy_dct.get(rpath_lvl_sp2))
                          if rpath_lvl_sp2 else None)

    # Set the ene thy info as a list of methods with coefficients
    ene_thy_info = []
    if isinstance(ene_lvl, str):
        ene_thy_info.append(
            [1.00, tinfo.from_dct(thy_dct.get(ene_lvl))])
    else:
        for lvl in ene_lvl:
            ene_thy_info.append(
                [lvl[0], tinfo.from_dct(thy_dct.get(thy_dct))])

    # Combine levels into a list
    es_levels = {
        'geo': (geo_lvl, geo_thy_info),
        'ene': (ene_lvl, ene_thy_info),
        'harm': (harm_lvl, harm_thy_info),
        'vpt2': (vpt2_lvl, vpt2_thy_info),
        'sym': (sym_lvl, sym_thy_info),
        'etrans': (etrans_lvl, etrans_thy_info),
        'tors': ([tors_lvl_sp, tors_lvl_scn],
                 [tors_sp_thy_info, tors_scn_thy_info]),
        'rpath': ([rpath_lvl_sp, rpath_lvl_scn, rpath_lvl_sp2],
                  [rpath_sp_thy_info, rpath_scn_thy_info, rpath_sp2_thy_info])
    }

    return es_levels
