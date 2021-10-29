""" Parses the `models.dat` input file for MechDriver that specifies all
    of the parameters used to construct the partition functions and
    master equations for the calculation of rate constants and thermochemistry.

"""

import automol
import ioformat
from mechanalyzer.inf import thy as tinfo
from mechlib.amech_io.parser._keywrd import defaults_with_dcts


# DCTS
MODKIN_REQ_LST = ('pressures', 'rate_temps', 'thermo_temps')
MODKIN_VAL_DCT = {
    'pressures': ((tuple,), (), None),
    'rate_temps': ((tuple,), (), None),
    'thermo_temps': ((tuple,), (), None),
    'well_extension_pressure': ((int, float), (), 1.0),
    'well_extension_temp': ((int, float), (), 600.0),
    'temp_unit': ((str,), (), 'K'),
    'pressure_unit': ((str,), (), 'atm'),
    'rate_fit': {
        'fit_method': ((str,), ('arrhenius', 'chebyshev'), 'arrhenius'),
        'pdep_temps': (((tuple,),), (), (500, 100)),
        'pdep_tol': ((float,), (), 20.0),
        'pdep_pval': ((float,), (), 1.0),
        'pdep_low': ((float,), (), None),
        'pdep_high': ((float,), (), None),
        'arrfit_dbltol': ((float,), (), 15.0),
        'arrfit_dblcheck': ((str,), ('max', 'avg'), 'max'),
        'troefit_params': ((tuple,), (), ('ts1', 'ts2', 'ts3', 'alpha')),
        'chebfit_tdeg': (int, (), 6),
        'chebfit_pdeg': (int, (), 4),
        'chebfit_tol': ((float,), (), 20.0)
    },
    'thermo_fit': {
        'ref_scheme': ((str,), ('basic', 'cbh0'), 'basic'),
        'ref_enes': ((str,), ('ANL0',), 'ANL0')
    },
    'glob_etransfer': {
        'lj': (None, (), 'estimate'),
        'edown': (None, (), 'estimate'),
        'mass': ((tuple,), (), None)
    }
}

MODPF_REQ_LST = ('pressures', 'rate_temps', 'thermo_temps')
MODPF_VAL_DCT = {
    'ene': {
        'lvl1': ((tuple,), (), None),
        'lvl2': ((tuple,), (), None)
    },
    'rot': {
        'mod': ((str,), ('rigid', 'vpt2'), 'rigid'),
        'vpt2lvl': ((str,), (), None)
    },
    'vib': {
        'mod': ((str,), ('harm', 'vpt2', 'tau'), 'harm'),
        'geolvl': ((str,), (), None),
        'vpt2lvl': ((str,), (), None),
    },
    'tors': {
        'mod': (
            (str,),
            ('rigid', '1dhr', '1dhrf', '1dhrfa', 'mdhr', 'mdhrv', 'tau'),
            'rigid'),
        'enelvl': ((str,), (), None),
        'geolvl': ((str,), (), None),
    },
    'symm': {
        'mod': ((str,), ('none', 'sampling', '1dhr'), 'none'),
        'geolvl': ((str,), (), None),
    },
    'rpath': {
        'enelvl': ((str,), (), None),
        'geolvl': ((str,), (), None),
    },
    'ts': {
        'nobar': ((str,), ('pst', 'rpvtst', 'vrctst'), 'pst'),
        'sadpt': ((str,), ('fixed', 'pst', 'rpvtst', 'vrctst'), 'fixed'),
        'rwells': ((str,), ('fake', 'find', 'none'), 'fake'),
        'pwells': ((str,), ('fake', 'find', 'none'), 'fake'),
        'tunnel': ((str,), ('none', 'eckart', 'sct'), 'eckart'),
        'etrans': ((str,), ('none', 'estimate', 'read'), 'estimate')
    }
}


# Build Basic Objects
def models_dictionary(mod_str, thy_dct):
    """ Parse the models.dat file string and construct the corresponding
        model dictionaries that are used by other parts of the code.

        :param mod_str: models.dat input file string
        :type mod_str: str
        :param thy_dct: information
        :type: dict[str:dict[str:str/float]]
    """

    # Format the models input to the kin and spc model dcts
    kin_blocks = ioformat.ptt.named_end_blocks(
        mod_str, 'kin', footer='kin')
    spc_blocks = ioformat.ptt.named_end_blocks(
        mod_str, 'spc', footer='spc')

    # Add defaults, check key-vals, and format each model dicts
    if kin_blocks is not None:
        kin_mod_dct = automol.util.dict_.merge_subdct(
            ioformat.ptt.keyword_dcts_from_blocks(kin_blocks),
            keep_subdct=True)
        for mod, dct in kin_mod_dct.items():
            if dct:  # if statement for empty global dcts from above fxn
                kin_mod_dct[mod] = _kin_model_build(dct)
    else:
        kin_mod_dct = None

    if spc_blocks is not None:
        spc_mod_dct = automol.util.dict_.merge_subdct(
            ioformat.ptt.keyword_dcts_from_blocks(spc_blocks),
            keep_subdct=True)
        for mod, dct in spc_mod_dct.items():
            if dct:  # if statement for empty global dcts from above fxn
                spc_mod_dct[mod] = _spc_model_build(dct, thy_dct)
    else:
        spc_mod_dct = None

    return kin_mod_dct, spc_mod_dct


# Convert objects
def _spc_model_build(spc_model_dct_i, thy_dct):
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

    # Have to check the dictionary to see if levels
    # in mod are in thy dct or it breaks

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
                thy_info = _format_lvl(val2)
                _new_dct[key2] = (val2, thy_info)
            else:
                _new_dct[key2] = val2

        new_dct2[key1] = _new_dct

    return new_dct2


def _kin_model_build(kin_mod_dct_i):
    """ set kin
    """

    # Set defaults
    new_kin_dct = automol.util.dict_.right_update(
        defaults_with_dcts(MODKIN_VAL_DCT), kin_mod_dct_i)

    # Check for correct input

    # Repartition ratefit key word into dcts for `ratefit` code
    old_ratefit = new_kin_dct['rate_fit']
    new_ratefit = {}
    for key, val in old_ratefit.items():
        if not any(x in key for x in ('pdep', 'arr', 'cheb', 'troe')):
            if val == 'arrhenius':  # keeps legacy working with the new form
                val = 'arr'
            elif val == 'chebyshev':
                val = 'cheb'
            new_ratefit[key] = val

    new_ratefit.update({
        'pdep_fit': {key.replace('pdep_', ''): val
                     for key, val in old_ratefit.items()
                     if 'pdep' in key},
        'arrfit_fit': {key.replace('arrfit_', ''): val
                       for key, val in old_ratefit.items()
                       if 'arrfit' in key},
        'chebfit_fit': {key.replace('chebfit_', ''): val
                        for key, val in old_ratefit.items()
                        if 'chebfit' in key},
        'troefit_fit': {key.replace('troefit_', ''): val
                        for key, val in old_ratefit.items()
                        if 'troefit' in key}
    })

    new_kin_dct.update({'rate_fit': new_ratefit})

    return new_kin_dct


# Handle model tasks
def extract_models(tsk):
    """ pull modesl from tasks
    """

    key_dct = tsk[-1]
    pes_mod = key_dct['kin_model']
    spc_mods = tuple(key_dct[key] for key in key_dct.keys()
                     if 'spc_mod' in key)

    return spc_mods, pes_mod


def split_model(mod):
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
    for char in mod:
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
    models.append(model)

    return (tuple(models), tuple(coeffs), tuple(operators))

# Functions to simplify accessing the model dictionaries
# def get_model(comp, spc_mod_dct_i):
#     """ get model for spc dct compoenent
#     """
#     return spc_mod_dct_i[comp]['mod']
# def thy_method_name(mod_lvl):
#     """ get the electronic strucure method
#     """
#     return mod_lvl[0]  # ?
# def thy_info_coeff(mod_lvl):
#     """ get the coefficient from the model
#     """
#     return mod_lvl[1][0]
# def thy_info_obj(mod_lvl):
#     """ get the name of the mod lvl
#     """
#     return mod_lvl[1][1]
