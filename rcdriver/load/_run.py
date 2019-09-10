from autoparse.find import all_captures
from autoparse.find import first_capture
from autoparse.pattern import capturing
from autoparse.pattern import zero_or_more
from autoparse.pattern import one_or_more
from autoparse.pattern import NONNEWLINE
from autoparse.pattern import NONSPACE
from autoparse.pattern import SPACE
from autoparse.pattern import WILDCARD
from autoparse.pattern import INTEGER
from autoparse.pattern import FLOAT

import logging
log   = logging.getLogger(__name__)


REQUIRED_KEYWORDS = [
    'reactants', 
    'jobs'
]

# Set the supported keywords for minima
SUPPORTED_KEYWORDS = [  
   'products',
   'paths',
   'references',
   'mc',
   'tsfind',
   'geom',
   'freq',
   'hr',
   'anharm',
   'var',
   'sp'
]
# FUNCTION TO READ IN A STRING FOR A SPECIFIC LEVEL # 

def get_task_str(inp_str, task):
    """ species input
    """
    # Set the species string pattern
    task_pattern = ('task' + one_or_more(SPACE) + task +  
                      capturing(one_or_more(WILDCARD, greedy=False)) +
                      'end')

    # Obtain the species string
    task_str = first_capture(task_pattern, inp_str)

    return task_str


# FUNCTIONS TO CHECK WHICH SPECIES AND KEYWORDS HAVE BEEN DEFINED

def get_defined_tasks(inp_str):
    """ gets a list which specifies what levels have been defined
    """
    levels_def_pattern =  ('task' + one_or_more(SPACE) + capturing(one_or_more(NONSPACE)))  

    defined_levels = all_captures(levels_def_pattern, inp_str)

    return defined_levels


def get_defined_keywords(inp_str, task):
    """ obtains the keywords defined in the input by the user
    """
    
    # Obtain the appropriate species string
    task_str = get_task_str(inp_str, task)

    # Get a list of the defined keywords
    defined_keys = []
    for line in task_str.splitlines():
        if '=' in line:
            tmp = line.strip().split('=')[0]
            defined_keys.append(tmp.strip())

    return defined_keys

def get_attr_call(key):
    calls = {'reactants' :  'get_key_array', 
             'products'  :  'get_key_array', 
             'paths'     :  'get_key_array', 
             'job'      :  'get_key_array', 
             'references' :  'get_key_array', 
              'mc'       :  'get_key_array',
              'tsfind'       :  'get_key_array',
              'geom'     :  'get_key_array',
              'freq'     :  'get_key_array',
              'hr'       :  'get_key_array',
              'anharm'   :  'get_key_array',
              'var'      :  'get_key_array',
              'sp'       :  'get_key_array'}
    if key in calls:
        return calls[key]
    else:
        log.warning('No method call for ooooo trait {}'.format(key))
        return 'error_call'

def error_call(s, task):
    log.warning('for task {}'.format(task))
    return None

def get_key_array(inp_str, task, key):
    # Find the value for a given key, make int
    pattern = (key + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                      capturing(one_or_more(NONNEWLINE)))
    # Obtain the appropriate species string
    patt_str = get_task_str(inp_str, task)
    # Obtain the charge string
    patt_str = first_capture(pattern, patt_str)
    assert patt_str is not None
    return patt_str.split()

