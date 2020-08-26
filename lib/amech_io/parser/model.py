""" Library of parser functions for the model file
"""

import sys
import copy
import autoparse.find as apf
from lib import filesys
from lib.amech_io.parser import ptt
from lib.amech_io.parser.keywords import MODEL_PF_SUPPORTED_DCT
from lib.amech_io.parser.keywords import MODEL_PF_DEFAULT_DCT


MODEL_INP = 'inp/models.dat'


# FUNCTION TO READ IN A STRING FOR A SPECIFIC LEVEL #

def read_models_sections(job_path):
    """ species input
    """
    mod_str = ptt.read_inp_str(job_path, MODEL_INP, remove_comments='#')

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
    pes_model_methods['global'] = glob_pes_model_methods
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
    spc_model_methods['global'] = glob_spc_model_methods
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
    therm_temps_str = apf.first_capture(
        ptt.paren_section('therm_temps'), model_str)
    rate_temps_str = apf.first_capture(
        ptt.paren_section('rate_temps'), model_str)
    pressures_str = apf.first_capture(
        ptt.paren_section('pressures'), model_str)
    etrans_str = apf.first_capture(ptt.paren_section('etransfer'), model_str)
    pdep_str = apf.first_capture(ptt.paren_section('pdep_fit'), model_str)
    tunit = apf.first_capture(ptt.keyword_pattern('tunit'), model_str)
    punit = apf.first_capture(ptt.keyword_pattern('punit'), model_str)
    fitm = apf.first_capture(ptt.keyword_pattern('fit_method'), model_str)
    ethr = apf.first_capture(
        ptt.keyword_pattern('dbl_arrfit_thresh'), model_str)
    echk = apf.first_capture(
        ptt.keyword_pattern('dbl_arrfit_check'), model_str)
    # assert temps_str is not None
    assert pressures_str is not None
    assert etrans_str is not None

    # Get the dictionary/values for each section and check them
    # Setting defaults
    therm_temps_dct, rate_temps_dct = {}, {}
    if temps_str is not None:
        therm_temps_dct = ptt.build_vals_lst(temps_str)
        rate_temps_dct = ptt.build_vals_lst(temps_str)
    if therm_temps_str is not None:
        therm_temps_dct = ptt.build_vals_lst(therm_temps_str)
    if rate_temps_str is not None:
        rate_temps_dct = ptt.build_vals_lst(rate_temps_str)

    pressures_dct = ptt.build_vals_lst(pressures_str)
    etransfer_dct = ptt.build_keyword_dct(etrans_str)
    pdep_dct = ptt.build_keyword_dct(pdep_str) if pdep_str is not None else {}
    tunit = ptt.set_value_type(tunit) if tunit is not None else 'K'
    punit = ptt.set_value_type(punit) if punit is not None else 'atm'
    fitm = ptt.set_value_type(fitm) if fitm is not None else 'arrhenius'
    ethr = ptt.set_value_type(ethr) if ethr is not None else 15.0
    echk = ptt.set_value_type(echk) if echk is not None else 'max'

    # Combine dcts into single model dct
    model_dct = {}
    model_dct['rate_temps'] = rate_temps_dct
    model_dct['therm_temps'] = therm_temps_dct
    model_dct['pressures'] = pressures_dct
    model_dct['etransfer'] = etransfer_dct
    model_dct['pdep_fit'] = pdep_dct
    model_dct['tunit'] = tunit
    model_dct['punit'] = punit
    model_dct['fit_method'] = fitm
    model_dct['dbl_arrfit_thresh'] = ethr
    model_dct['dbl_arrfit_check'] = echk
    # model_dct['nasa_temp_ranges'] = [200, 1000, 3000]

    # Check the dct

    return model_dct


def build_spc_model_keyword_dct(model_str):
    """ Build a dictionary for all the models keywords
    """
    # Grab the various sections required for each model
    pf_str = apf.first_capture(ptt.paren_section('pf'), model_str)
    es_str = apf.first_capture(ptt.paren_section('es'), model_str)
    vrctst_str = apf.first_capture(ptt.paren_section('vrctst'), model_str)
    opts_str = apf.first_capture(ptt.paren_section('options'), model_str)
    if pf_str is None:
        print('*ERROR: pf section is not defined')
        sys.exit()
    if es_str is None:
        print('*ERROR: es section is not defined')
        sys.exit()

    # Get the dictionary for each section and check them
    pf_dct = ptt.build_keyword_dct(pf_str)
    es_dct = ptt.build_keyword_dct(es_str)
    opts_dct = ptt.build_keyword_dct(opts_str)
    if 'ref_enes' not in opts_dct:
        opts_dct['ref_enes'] = 'ATcT'
    if 'ref_scheme' not in opts_dct:
        opts_dct['ref_scheme'] = 'basic'
    # if vrctst_str is not None:
    #     vrctst_dct = ptt.build_keyword_dct(vrctst_str)
    # else:
    #     vrctst_dct = _vrc_dct()

    # Combine dcts into single model dct
    model_dct = {}
    model_dct['pf'] = pf_dct
    model_dct['es'] = es_dct
    # model_dct['vrctst'] = vrctst_dct
    model_dct['options'] = opts_dct

    # Check the dct
    check_spc_model_dct(model_dct)

    return model_dct


# def check_pes_model_dct(model_dct):
#     """ Make sure the models dictionary keywords are all correct
#     """
#     pass


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


def set_default_pf(dct):
    """ set defaults
    """
    new_dct = copy.deepcopy(dct)
    for key, val in MODEL_PF_DEFAULT_DCT.items():
        if key not in new_dct:
            new_dct[key] = val

    return new_dct


def pf_model_info(pf_model):
    """ Set the PF model list based on the input
    """
    print('pfmodel', pf_model)
    rot_model = pf_model['rot'] if 'rot' in pf_model else 'rigid'
    tors_model = pf_model['tors'] if 'tors' in pf_model else 'rigid'
    vib_model = pf_model['vib'] if 'vib' in pf_model else 'harm'
    sym_model = pf_model['sym'] if 'sym' in pf_model else 'none'

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
        'rwells': rwells_model,
        'pwells': pwells_model
    }
    print('pfmodels', pf_models)

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

    # Torsions and rxn paths which needs a reference for itself
    tors_lvl_sp = es_model['tors'][0] if 'tors' in es_model else None
    tors_lvl_scn = es_model['tors'][1] if 'tors' in es_model else None
    rpath_lvl_sp = es_model['rpath'][0] if 'rpath' in es_model else None
    rpath_lvl_scn = es_model['rpath'][1] if 'rpath' in es_model else None

    # Set the theory info objects
    geo_thy_info = filesys.inf.get_thy_info(geo_lvl, thy_dct)
    harm_thy_info = filesys.inf.get_thy_info(harm_lvl, thy_dct)
    vpt2_thy_info = (filesys.inf.get_thy_info(vpt2_lvl, thy_dct)
                     if vpt2_lvl else None)
    sym_thy_info = (filesys.inf.get_thy_info(sym_lvl, thy_dct)
                    if sym_lvl else None)
    tors_sp_thy_info = (filesys.inf.get_thy_info(tors_lvl_sp, thy_dct)
                        if tors_lvl_sp else None)
    tors_scn_thy_info = (filesys.inf.get_thy_info(tors_lvl_scn, thy_dct)
                         if tors_lvl_scn else None)
    rpath_sp_thy_info = (filesys.inf.get_thy_info(rpath_lvl_sp, thy_dct)
                         if rpath_lvl_sp else None)
    rpath_scn_thy_info = (filesys.inf.get_thy_info(rpath_lvl_scn, thy_dct)
                          if rpath_lvl_scn else None)

    # Set the ene thy info as a list of methods with coefficients
    ene_thy_info = []
    if isinstance(ene_lvl, str):
        ene_thy_info.append(
            [1.00, filesys.inf.get_thy_info(ene_lvl, thy_dct)])
    else:
        for lvl in ene_lvl:
            ene_thy_info.append(
                [lvl[0], filesys.inf.get_thy_info(lvl[1], thy_dct)])

    # Combine levels into a list
    es_levels = {
        'geo': (geo_lvl, geo_thy_info),
        'ene': (ene_lvl, ene_thy_info),
        'harm': (harm_lvl, harm_thy_info),
        'vpt2': (vpt2_lvl, vpt2_thy_info),
        'sym': (sym_lvl, sym_thy_info),
        'tors': ([tors_lvl_sp, tors_lvl_scn],
                 [tors_sp_thy_info, tors_scn_thy_info]),
        'rpath': ([rpath_lvl_sp, rpath_lvl_scn],
                  [rpath_sp_thy_info, rpath_scn_thy_info])
    }

    return es_levels


def _vrc_dct():
    """ Build VRC dict
    """
    return {
        'fortran_compiler': 'gfortran',
        'spc_name': 'mol',
        'memory': 4.0,
        'r1dists_lr': [8., 6., 5., 4.5, 4.],
        'r1dists_sr': [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2],
        'r2dists_sr': [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2],
        'd1dists': [0.01, 0.5, 1.],
        'd2dists': [0.01, 0.5, 1.],
        'conditions': {},
        'nsamp_max': 2000,
        'nsamp_min': 50,
        'flux_err': 10,
        'pes_size': 2,
        'exe_path': '/blues/gpfs/home/sjklipp/bin/molpro'
    }
