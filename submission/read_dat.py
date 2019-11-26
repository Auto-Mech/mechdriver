""" read the params from a dat file
"""

from types import SimpleNamespace

# setting SORT_RXNS to False leads to missing channels
# for now just leave them sorted

# Set each keyword to their default values
PARAMS = {
    'DATA_PATH': '/home/sjklipp/PACC/mech_test',
    'HIND_INC': 30.,
    'MC_NSAMP0': [True, 6, 1, 3, 100],
    'REF_MOLS': 'basic',
    'ENE_COEFF': [1.],
    'PESNUMS': 'all',
    'CHANNELS': 'all',
    'CHECK_STEREO': False,
    'SORT_RXNS': True,
    'RAD_RAD_SORT': True,
    'RAD_RAD_TS': 'vtst',
    'PST_PARAMS': [1.0, 6],
    'RUN_THERMO': False,
    'RUN_RATES': True,
    'OPT_ES': True,
    'OPT_MESS': False,
    'OPT_THERMO': False,
    'OPT_ALLPF': False,
    'RUN_ES': True,
    'RUN_SPECIES': True,
    'RUN_MESS': True,
    'RUN_RATES_OPT': True,
    'OPT_LVL0': 'lvl_wbs',
    'OPT_LVL1': 'lvl_wbm',
    'OPT_LVL2': 'lvl_b2t',
    'SCAN_LVL1': 'lvl_wbs',
    'SP_LVL1': 'cc_lvl_df',
    'SP_LVL2': 'cc_lvl_tf',
    'SP_LVL3': 'cc_lvl_qf',
    'MULTI_LVL': 'mlvl_cas_dz',
    'TSK_INFO_LST': [
        ['find_geom', 'lvl_wbs', 'lvl_wbs', False],
        ['conf_samp', 'lvl_wbs', 'lvl_wbs', False],
        ['conf_hess', 'lvl_wbs', 'lvl_wbs', False],
        ['sym_samp', 'lvl_wbs', 'lvl_wbs', False],
        ['hr_scan', 'lvl_wbs', 'lvl_wbs', False],
        ['find_ts', 'lvl_wbs', 'lvl_wbs', False],
        ['conf_samp', 'lvl_wbs', 'lvl_wbs', False],
        ['conf_hess', 'lvl_wbs', 'lvl_wbs', False],
        ['sym_samp', 'lvl_wbs', 'lvl_wbs', False],
        ['hr_scan', 'lvl_wbs', 'lvl_wbs', False],
        ['conf_energy', 'cc_lvl_df', 'lvl_wbs', False],
        ['sym_samp', 'lvl_wbs', 'lvl_wbs', False],
        ['conf_vpt2', 'lvl_wbm', 'lvl_wbm', False],
    ],
    'OVERWRITE': False,
    'EXP_FACTOR': 150.0,
    'EXP_POWER': 0.85,
    'EXP_CUTOFF': 15.,
    'EPS1': 100.0,
    'EPS2': 200.0,
    'SIG1': 6.,
    'SIG2': 6.,
    'MASS1': 15.0,
    'RUN_PREFIX': '/lcrc/project/PACC/run',
    'SAVE_PREFIX': '/lcrc/project/PACC/save'
}


def params(input_file):
    """ read the run parameters from a file
    """

    # Read the parameters from the input file ignoring task list keywords
    with open(input_file, 'r') as inp_file:
        for line in inp_file:
            if '#' not in line and line.strip() != '':
                param_vals = line.strip().split('=')
                if param_vals:
                    keyword, value = _format_param_vals(param_vals)
                    if 'TSK_' not in keyword:
                        PARAMS[keyword] = value

    # Build a name space object for the parameters
    set_params = SimpleNamespace(**PARAMS)

    # Set helpful options lists using the values of the parameters
    setattr(set_params, 'OPTIONS_THERMO',
            [set_params.OPT_ES, set_params.OPT_MESS,
             set_params.OPT_THERMO, set_params.OPT_ALLPF])
    setattr(set_params, 'OPTIONS_RATE',
            [set_params.RUN_ES, set_params.RUN_SPECIES,
             set_params.RUN_MESS, set_params.RUN_RATES_OPT])

    # Now read the parameters used to set the TSK_INFO_LST
    tsk_info_lst = []
    with open(input_file, 'r') as inp_file:
        for line in inp_file:
            if '#' not in line and line.strip() != '':
                param_vals = line.strip().split('=')
                if param_vals:
                    keyword, value_lst = _format_param_vals(param_vals)
                    if 'TSK_' in keyword:
                        tsk_info_lst.append(
                            _format_task_lst(set_params, keyword, value_lst))

    # Set the task list as an attribute of the set_params namespace
    setattr(set_params, 'TSK_INFO_LST', tsk_info_lst)

    return set_params


def _format_param_vals(pvals):
    """ format param vals string
    """
    [keyword, value] = pvals

    frmtd_keyword = keyword.strip().upper()

    value = value.strip()
    # Format values if it is a list (of string(s), boolean(s), int(s))
    if all(sym in value for sym in ('[', ']', ',')):
        value = value.replace('[', '').replace(']', '')
        value = value.split(',')
        frmtd_value = []
        # Set string in list to boolean or integer if needed
        for elm in value:
            frmtd_value.append(_set_value_type(elm.strip()))
    else:
        # Format values if it has singular value
        frmtd_value = _set_value_type(value)

    return frmtd_keyword, frmtd_value


def _format_task_lst(set_params, keyword, value_lst):
    """ format keyword = [vals] string to the tsk_lst keywords
    """
    frmtd_keyword = keyword.split('TSK_')[1].lower()

    setvals = [getattr(set_params, value.strip())
               for value in value_lst]
    tsk_lst = [frmtd_keyword, setvals[0], setvals[1], setvals[2]]

    return tsk_lst


def _set_value_type(value):
    """ set type of value
        right now we handle True/False boolean, int, float, and string
    """

    if value == 'True':
        frmtd_value = True
    elif value == 'False':
        frmtd_value = False
    elif value.isdigit():
        frmtd_value = int(value)
    elif '.' in value and value.replace('.', '').replace('-','').isdigit():
        frmtd_value = float(value)
    else:
        frmtd_value = value

    return frmtd_value
