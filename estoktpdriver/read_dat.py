""" read the params from a dat file
"""

from types import SimpleNamespace

# Set each keyword to their default values
PARAMS = {
    'SORT_RXN': True,
    'HIND_INC': 30.,
    'RUN_THERMO': False,
    'RUN_RATES': True,
    'CHECK_STEREO': False,
    'MC_NSAMP0': [True, 6, 1, 3, 100],
    'OPT_ES': True,
    'OPT_MESS': False,
    'OPT_THERMO': False,
    'OPT_ALLPF': False,
    'REF_MOLS': 'basic',
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
    'OVERWRITE': False,
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
    ]
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
            elm = elm.strip()
            if elm == 'true':
                elm = True
            elif elm == 'false':
                elm = False
            elif elm.isdigit():
                elm = int(elm)
            elif '.' in value and elm.replace('.', '').isdigit():
                elm = float(elm)
            else:
                elm = elm.replace("'", '')

            frmtd_value.append(elm)
    # Format values if it is just a boolean
    elif value == 'true':
        frmtd_value = True
    elif value == 'false':
        frmtd_value = False
    # Simply set if it is a string
    else:
        frmtd_value = value
        frmtd_value = value.replace("'", '')

    return frmtd_keyword, frmtd_value


def _format_task_lst(set_params, keyword, value_lst):
    """ format keyword = [vals] string to the tsk_lst keywords
    """
    frmtd_keyword = keyword.split('TSK_')[1].lower()

    setvals = [getattr(set_params, value.strip())
               for value in value_lst]
    tsk_lst = [frmtd_keyword, setvals[0], setvals[1], setvals[2]]

    return tsk_lst
