"""
general task prints
"""

from mechlib.amech_io.printer._print import message
from mechlib.amech_io.printer._lib import obj


TASK_STR = """Task
  {}

Species
  Name: {}
  SMILES: {}
"""


def task_header(tsk, spc_name):
    """ a
    """
    obj('vspace')
    obj('line_dash')
    message('Task:', tsk, spc_name, newline=1)


def task_footer():
    """ a
    """
    message(
        '\n'+
        '--- END TASK ----------------------------------------' +
        '--------------------------------------'
    )


def keyword_list(es_keyword_dct, thy_dct=None):
    """ a
    """
    message('Options for electronic structure task:', newline=1)
    for key, val in es_keyword_dct.items():
        method_str = ''
        if key in ('inplvl', 'runlvl'):
            method_dct = thy_dct.get(es_keyword_dct[key])
            method_str = f'({method_dct["method"]}/{method_dct["basis"]})'
        message(f'{key}: {val}    ' + method_str)
    obj('vspace')


def output_task_header(tsk):
    """ a
    """
    obj('vspace')
    obj('line_dash')
    message('Print Property:', tsk, newline=1)


def output_keyword_list(es_keyword_dct, thy_dct=None):
    """ a
    """
    message('Electronic structure level for property:', newline=1)
    for key, val in es_keyword_dct.items():
        method_str = ''
        if key in ('inplvl', 'runlvl'):
            method_dct = thy_dct.get(es_keyword_dct[key])
            method_str = f'({method_dct["method"]}/{method_dct["basis"]})'
        message(f'{key}: {val}    ' + method_str)
    obj('vspace')


def messpf(statement, path=None):
    """ a
    """
    # obj('vspace')
    if statement == 'write_header':
        obj('line_dash')
        message('Preparing MESSPF input files for all species', newline=1)
    elif statement == 'input_string':
        message('MESSPF Input String:')
        obj('vspace')
        obj('vspace')
    elif statement == 'run_header':
        obj('line_dash')
        message('Running MESSPF calculations for all species', newline=1)
    elif statement == 'write_file':
        message('Writing MESS input file...')
        message(f' - Path: {path}')
    elif statement == 'write_output':
        message('Writing MESS Output file...')
        message(f' - Path: {path}')
    elif statement == 'run_file':
        message('Running MESS input file...')
        message(f' - Path: {path}')
    elif statement == 'global_header':
        message('Preparing global keywords section for MESS input...')
        message(' - Using temperatures and pressures defined by user')
        message(' - Using internal AutoMech defaults for other MESS keywords:')
    elif statement == 'transfer_section':
        message('Preparing energy transfer section for MESS input...')
    elif statement == 'well_section':
        message('- Determining reference well species...')
    elif statement == 'bath_section':
        message('- Determining information for the bath species...')
    elif statement == 'channel_section':
        message('Preparing reaction channel section for MESS input... ')


def nasa(statement, spc_name=None, path=None):
    """ a
    """
    obj('vspace')
    if statement == 'header':
        obj('line_dash')
        message(
            'Running Thermochemistry calculations for all species', newline=1)
    elif statement == 'calculate':
        message(
            f'Starting NASA polynomials calculation for {spc_name}')
    elif statement == 'fit':
        message(
            'Attempting to fit NASA polynomials from',
            '200-1000 and 1000-3000 K ranges using\n',
            f'temps from MESSPF file:\n {path}.')
