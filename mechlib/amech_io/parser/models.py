""" Library of parser functions for the model file

    I think it this is very out of wack with current set-up
"""

import automol
import ioformat
from mechanalyzer.inf import thy as tinfo
# from mechlib.amech_io.parser._keywrd import defaults_from_val_dct
from mechlib.amech_io.parser._keywrd import defaults_with_dcts
from mechlib.amech_io.parser._keywrd import MODPF_VAL_DCT


# Build Basic Objects
def models_dictionary(mod_str):
    """ Parse the models.dat file
    """

    # Format the models input to the kin and spc model dcts
    kin_blocks = ioformat.ptt.named_end_blocks(mod_str, 'kin')
    spc_blocks = ioformat.ptt.named_end_blocks(mod_str, 'spc')

    kin_mod_dct = automol.util.dict_.merge_subdct(
        ioformat.ptt.keyword_dcts_from_blocks(kin_blocks), keep_subdct=True)
    spc_mod_dct = automol.util.dict_.merge_subdct(
        ioformat.ptt.keyword_dcts_from_blocks(spc_blocks), keep_subdct=True)

    # Assess if the model.dat input is valid
    # kin_mod_dct = automol.util.dict_.right_update(
    #     some_defaults_fxn(MODKIN_VAL_DCT), kin_mod_dct)

    return kin_mod_dct, spc_mod_dct


# Convert objects
def pf_model_info(spc_model_dct):
    """ Set the PF model list based on the input
        Combine with es

        {'vib': {'model': '', 'geolvl_info': (), ...}}
    """

    if 'wells' in spc_model_dct['ts']:
        spc_model_dct['ts']['rwells'] = spc_model_dct['ts']['wells']
        spc_model_dct['ts']['pwells'] = spc_model_dct['ts']['wells']
        spc_model_dct['ts'].pop('wells')

    new_dct = automol.util.dict_.right_update(
      defaults_with_dcts(MODPF_VAL_DCT), spc_model_dct)

    return new_dct


def pf_level_info(spc_model_dct, thy_dct):
    """ Set the model info

        currently breaks for composite methods, e.g. energies
        maybe just have a function grab all of them and calculate it with the
        levels later
    """

    def _format_val(lvl_val):
        """ format weird energy calls
        """
        if isinstance(lvl_val, str):
            val_inf = [1.00, tinfo.from_dct(thy_dct.get(ene_lvl))]
        else:
            val_inf = [lvl[0], tinfo.from_dct(thy_dct.get(thy_dct))]

        return val_inf


    new_dct = {}
    for key1, val1 in spc_model_dct.items():
        _new_dct = {}
        for key2, val2 in val1.items():
            if 'lvl' in key2:
                thy_info = _format_lvl(val2)
                # if isinstance(val2, str):
                #     thy_info = (tinfo.from_dct(thy_dct.get(val2))
                #                 if val2 else None)
                # else:
                #     pass
                #     # do the theory list setting with coefficients
                _new_dct[key2] = thy_info
            else:
                _new_dct[key2] = val2

        new_dct[key1] = _new_dct


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
