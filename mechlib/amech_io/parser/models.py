""" Library of parser functions for the model file

    I think it this is very out of wack with current set-up
"""

import automol
import ioformat
from mechanalyzer.inf import thy as tinfo
from mechlib.amech_io.parser._keywrd import defaults_with_dcts


# DCTS
MODKIN_VAL_DCT = {
    'pressures': (tuple, (), None),
    'rate_temps': (tuple, (), None),
    'thermo_temps': (tuple, (), None),
    'rate_fit': {
        'fit_method': (str, ('arrhenius', 'chebyshev'), 'arrhenius'),
        'pdep_temps': (tuple, (), (500, 100)),
        'pdep_tol': (float, (), 20.0),
        'pdep_pval': (float, (), 1.0),
        'pdep_low': (float, (), None),
        'pdep_high': (float, (), None),
        'arr_dbl_tol': (float, (), 15.0),
        'troe_param_fit_list': (
            tuple, (), ('ts1', 'ts2', 'ts3', 'alpha'))
    },
    'thermo_fit': {
        'ref_scheme': (str, ('basic', 'cbh0'), 'basic'),
        'ref_enes': (str, ('ANL0',), 'ANL0')
    },
    'glob_etransfer': {
        'lj': (None, (), 'estimate'),
        'edown': (None, (), 'estimate'),
        'mass': (tuple, (), None)
    }
}

MODPF_VAL_DCT = {
    'ene': {
        'lvl1': (tuple, (), None),
        'lvl2': (tuple, (), None)
    },
    'rot': {
        'mod': (str, ('rigid', 'vpt2'), 'rigid'),
        'vpt2lvl': (str, (), None)
    },
    'vib': {
        'mod': (str, ('harm', 'vpt2', 'tau'), 'harm'),
        'geolvl': (str, (), None),
        'vpt2lvl': (str, (), None),
    },
    'tors': {
        'mod': (
            str, ('rigid', '1dhr', '1dhrf', '1dhrfa', 'mdhr', 'mdhrv', 'tau'),
            'rigid'),
        'enelvl': (str, (), None),
        'geolvl': (str, (), None),
    },
    'symm': {
        'mod': (str, ('none', 'sampling', '1dhr'), 'none'),
        'geolvl': (str, (), None),
    },
    'rpath': {
        'enelvl': (str, (), None),
        'geolvl': (str, (), None),
    },
    'ts': {
        'nobar': (str, ('pst', 'rpvtst', 'vrctst'), 'pst'),
        'sadpt': (str, ('fixed', 'pst', 'rpvtst', 'vrctst'), 'fixed'),
        'rwells': (str, ('fake', 'find', 'none'), 'fake'),
        'pwells': (str, ('fake', 'find', 'none'), 'fake'),
        'tunnel': (str, ('none', 'eckart', 'sct'), 'eckart'),
        'etrans': (str, ('none', 'estimate', 'read'), 'estimate')
    }
}


# Build Basic Objects
def models_dictionary(mod_str, thy_dct):
    """ Parse the models.dat file
    """

    # Format the models input to the kin and spc model dcts
    kin_blocks = ioformat.ptt.named_end_blocks(mod_str, 'kin')
    spc_blocks = ioformat.ptt.named_end_blocks(mod_str, 'spc')

    kin_mod_dct = automol.util.dict_.merge_subdct(
        ioformat.ptt.keyword_dcts_from_blocks(kin_blocks), keep_subdct=True)
    spc_mod_dct = automol.util.dict_.merge_subdct(
        ioformat.ptt.keyword_dcts_from_blocks(spc_blocks), keep_subdct=True)

    print('kin', kin_mod_dct)
    print('spc', spc_mod_dct)

    # Add defaults and format each kin dictionary
    for mod, dct in kin_mod_dct.items():
        kin_mod_dct[mod] = automol.util.dict_.right_update(
            defaults_with_dcts(MODKIN_VAL_DCT), kin_mod_dct[mod])

    # Add defaults and format each spc dictionary
    for mod, dct in spc_mod_dct.items():
        spc_mod_dct[mod] = _spc_model_defaults(dct, thy_dct)

    return kin_mod_dct, spc_mod_dct


# Convert objects
def _spc_model_defaults(spc_model_dct_i, thy_dct):
    """ Set the PF model list based on the input
        Combine with es

        {'vib': {'model': '', 'geolvl_info': (), ...}}

        Take all of the spc model dcts or just one?
    """

    # Expand the wells keyword if it has been set
    if 'wells' in spc_model_dct_i['ts']:
        spc_model_dct_i['ts']['rwells'] = spc_model_dct_i['ts']['wells']
        spc_model_dct_i['ts']['pwells'] = spc_model_dct_i['ts']['wells']
        spc_model_dct_i['ts'].pop('wells')

    # Set the defaults for the dct
    new_dct = automol.util.dict_.right_update(
      defaults_with_dcts(MODPF_VAL_DCT), spc_model_dct_i)

    # Have to check the dictionary to see if levels in mod are in thy dct or it breaks

    # Format the keys that are thy levels into info objects for later use
    def _format_lvl(lvl_val):
        """ format weird energy calls
        """
        if isinstance(lvl_val, str):
            val_inf = (1.00, tinfo.from_dct(thy_dct.get(lvl_val)))
        else:
            val_inf = (lvl_val[0], tinfo.from_dct(thy_dct.get(lvl_val)))

        return val_inf

    new_dct2 = {}
    for key1, val1 in new_dct.items():
        _new_dct = {}
        for key2, val2 in val1.items():
            if 'lvl' in key2 and val2 is not None:
                print('key, val', key2, val2)
                thy_info = _format_lvl(val2)
                _new_dct[key2] = (val2, thy_info)
            else:
                _new_dct[key2] = val2

        new_dct2[key1] = _new_dct

    return new_dct2


def mult_models(mods, spc_model_dct, thy_dct):
    """ Build dictionaries with models
    """

    pf_levels = {}
    pf_models = {}
    for mod in mods:
        pf_levels[mod] = pf_level_info(spc_model_dct[mod]['es'], thy_dct)
        pf_models[mod] = pf_model_info(spc_model_dct[mod]['pf'])
        pf_models[mod]['ref_scheme'] = (
            spc_model_dct[mod]['options']['ref_scheme']
            if 'ref_scheme' in spc_model_dct[mod]['options'] else 'none')
        pf_models[mod]['ref_enes'] = (
            spc_model_dct[mod]['options']['ref_enes']
            if 'ref_enes' in spc_model_dct[mod]['options'] else 'none')

    return pf_levels, pf_models


def split_model(model):
    """ Take a model given by a set of operations and split it into a set of models
        and operators

        model = 3*pf1-pf2 ->
            model = (models, coefs, operators)
            model = ((pf1, pf2), (3, 1), (multiple, subtract))
    """

    # Dictionary to map symbols to strings used later in the code
    op_dct = {'*': 'multiply', '+': 'add', '/': 'divide', '-': 'substract'}

    # Break the down the model string into
    # constituent models, coefficients, and operators
    coeffs, operators, models = [], [], []
    coeff, model = '', ''
    for char in model:
        if char == '.' or char.isdigit():
            coeff += char
        elif char.isalpha():
            model += char
        elif char in op_dct:
            operators.append(op_dct[char])
            if coeff:
                coeffs.append(float(coeff))
            else:
                coeffs.append(1.0)
            models.append(model)
            coeff = ''
            model = ''
    if coeff:
        coeffs.append(float(coeff))
    else:
        coeffs.append(1.0)

    return (tuple(models), tuple(coeffs), tuple(operators))
