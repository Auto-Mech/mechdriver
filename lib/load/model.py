""" Library of reader functions for the model file
"""

import copy
import autoparse.find as apf
from lib.load import ptt
from lib.filesystem import inf as finf
from lib.load.keywords import MODEL_PF_SUPPORTED_DCT
from lib.load.keywords import MODEL_PF_DEFAULT_DCT


MODEL_INP = 'inp/models.dat'


# FUNCTION TO READ IN A STRING FOR A SPECIFIC LEVEL #

def read_models_sections(job_path):
    """ species input
    """
    mod_str = ptt.read_inp_str(job_path, MODEL_INP)

    # Build a dictionary for the PES models
    pes_model_sections = apf.all_captures(
        ptt.end_section_wname2('pes_model'), mod_str)
    assert pes_model_sections is not None

    pes_model_methods = {}
    glob_pes_model_methods = {}
    for section in pes_model_sections:
        if section[0] == 'global':
            name = section[0]
            keyword_dct = build_pes_model_keyword_dct(section[1])
            glob_pes_model_methods = keyword_dct
            break
    for section in pes_model_sections:
        if section[0] != 'global':
            name = section[0]
            keyword_dct = build_pes_model_keyword_dct(section[1])
            pes_model_methods[name] = combine_glob_and_spc_dct(
                glob_pes_model_methods, keyword_dct)

    # Build a dictionary for the spc models
    spc_model_sections = apf.all_captures(
        ptt.end_section_wname2('spc_model'), mod_str)
    assert spc_model_sections is not None

    spc_model_methods = {}
    glob_spc_model_methods = {}
    for section in spc_model_sections:
        if section[0] == 'global':
            name = section[0]
            keyword_dct = build_spc_model_keyword_dct(section[1])
            glob_spc_model_methods = keyword_dct
            break
    for section in spc_model_sections:
        if section[0] != 'global':
            name = section[0]
            keyword_dct = build_spc_model_keyword_dct(section[1])
            spc_model_methods[name] = combine_glob_and_spc_dct(
                glob_spc_model_methods, keyword_dct)

    return pes_model_methods, spc_model_methods


def combine_glob_and_spc_dct(glob_dct, spc_dct):
    """ stuff
    """
    new_dct = glob_dct.copy()
    for skey, sval in spc_dct.items():
        if skey in glob_dct:
            new_dct[skey] = sval

    return new_dct


def build_pes_model_keyword_dct(model_str):
    """ Build a dictionary for all the models keywords
    """
    # Grab the various sections required for each model
    temps_str = apf.first_capture(ptt.paren_section('temps'), model_str)
    pressures_str = apf.first_capture(
        ptt.paren_section('pressures'), model_str)
    etrans_str = apf.first_capture(ptt.paren_section('etransfer'), model_str)
    assert temps_str is not None
    assert pressures_str is not None
    assert etrans_str is not None

    # Get the dictionary for each section and check them
    temps_dct = ptt.build_vals_lst(temps_str)
    pressures_dct = ptt.build_vals_lst(pressures_str)
    etransfer_dct = ptt.build_keyword_dct(etrans_str)

    # Combine dcts into single model dct
    model_dct = {}
    model_dct['temps'] = temps_dct
    model_dct['pressures'] = pressures_dct
    model_dct['etransfer'] = etransfer_dct

    return model_dct


def build_spc_model_keyword_dct(model_str):
    """ Build a dictionary for all the models keywords
    """
    # Grab the various sections required for each model
    pf_str = apf.first_capture(ptt.paren_section('pf'), model_str)
    es_str = apf.first_capture(ptt.paren_section('es'), model_str)
    options_str = apf.first_capture(ptt.paren_section('options'), model_str)
    vrctst_str = apf.first_capture(ptt.paren_section('vrctst'), model_str)
    assert pf_str is not None
    assert es_str is not None
    assert options_str is not None

    # Get the dictionary for each section and check them
    pf_dct = ptt.build_keyword_dct(pf_str)
    es_dct = ptt.build_keyword_dct(es_str)
    options_dct = ptt.build_keyword_dct(options_str)
    if vrctst_str is not None:
        vrctst_dct = ptt.build_keyword_dct(vrctst_str)

    # Check the PF dct and set defaults
    assert check_model_dct(pf_dct)
    pf_dct = set_default_pf(pf_dct)

    # Combine dcts into single model dct
    model_dct = {}
    model_dct['pf'] = pf_dct
    model_dct['es'] = es_dct
    model_dct['options'] = options_dct
    if vrctst_str is not None:
        model_dct['vrctst'] = vrctst_dct

    return model_dct


def check_model_dct(dct):
    """ Make sure the models dictionary keywords are all correct
    """
    proper_def = True
    for key, val in dct.items():
        if key in MODEL_PF_SUPPORTED_DCT:
            if val not in MODEL_PF_SUPPORTED_DCT[key]:
                proper_def = False
        else:
            proper_def = False

    return proper_def


def set_default_pf(dct):
    """ set defaults
    """
    new_dct = copy.deepcopy(dct)
    for key, val in MODEL_PF_DEFAULT_DCT.items():
        if key not in new_dct:
            new_dct[key] = val

    return new_dct


def set_pf_model_info(pf_model):
    """ Set the PF model list based on the input
    """
    tors_model = pf_model['tors'] if pf_model['tors'] else 'RIGID'
    vib_model = pf_model['vib'] if pf_model['vib'] else 'HARM'
    sym_model = pf_model['sym'] if pf_model['sym'] else ''

    pf_models = [tors_model, vib_model, sym_model]

    return pf_models


def set_es_model_info(es_model, thy_dct):
    """ Set the model info
    """
    # Read the ES models from the model dictionary
    geo_lvl = es_model['geo'] if es_model['geo'] else None
    ene_lvl = es_model['ene'] if es_model['ene'] else None
    harm_lvl = es_model['harm'] if es_model['harm'] else None
    vpt2_lvl = es_model['vpt2'] if es_model['vpt2'] else None
    sym_lvl = es_model['sym'] if es_model['sym'] else None

    # Torsional Scan which needs a reference for itself
    tors_lvl_sp = es_model['tors'][0] if es_model['tors'] else None
    tors_lvl_scn = es_model['tors'][1] if es_model['tors'] else None

    # Set the theory info objects
    geo_thy_info = finf.get_thy_info(geo_lvl, thy_dct)
    ene_thy_info = finf.get_thy_info(ene_lvl, thy_dct)
    harm_thy_info = finf.get_thy_info(harm_lvl, thy_dct)
    vpt2_thy_info = (finf.get_thy_info(vpt2_lvl, thy_dct)
                     if vpt2_lvl else None)
    sym_thy_info = (finf.get_thy_info(sym_lvl, thy_dct)
                    if sym_lvl else None)
    tors_sp_thy_info = (finf.get_thy_info(tors_lvl_sp, thy_dct)
                        if tors_lvl_sp else None)
    tors_scn_thy_info = (finf.get_thy_info(tors_lvl_scn, thy_dct)
                         if tors_lvl_scn else None)

    # Set the ene thy info as a list of methods with coefficients
    ene_thy_info = ()
    if isinstance(ene_lvl, str):
        ene_thy_info += (1.00, finf.get_thy_info(ene_lvl, thy_dct))
    else:
        for lvl in ene_lvl:
            ene_thy_info += (lvl[0], finf.get_thy_info(lvl[1], thy_dct))

    # Combine levels into a list
    es_levels = [
        geo_thy_info,
        ene_thy_info,
        harm_thy_info,
        vpt2_thy_info,
        sym_thy_info,
        [tors_sp_thy_info, tors_scn_thy_info]
    ]

    return es_levels
