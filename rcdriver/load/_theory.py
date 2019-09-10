""" library of reader functions for the theory file
"""

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
    'program', 
    'method', 
    'basis'
]

# Set the supported keywords for minima
SUPPORTED_KEYWORDS = [  
   'ncycles',
   'orb_res',
   'mem'    ,
   'nprocs' ,
   'econv'  ,
   'gconv'  
]
# FUNCTION TO READ IN A STRING FOR A SPECIFIC LEVEL # 

def get_lvl_str(inp_str, lvl):
    """ species input
    """
    # Set the species string pattern
    lvl_pattern = ('level' + one_or_more(SPACE) + lvl +  
                      capturing(one_or_more(WILDCARD, greedy=False)) +
                      'end')

    # Obtain the species string
    lvl_str = first_capture(lvl_pattern, inp_str)

    return lvl_str


# FUNCTIONS TO CHECK WHICH SPECIES AND KEYWORDS HAVE BEEN DEFINED

def get_defined_lvls(inp_str):
    """ gets a list which specifies what levels have been defined
    """
    levels_def_pattern =  ('level' + one_or_more(SPACE) + capturing(one_or_more(NONSPACE)))  

    defined_levels = all_captures(levels_def_pattern, inp_str)

    return defined_levels


def get_defined_keywords(inp_str, lvl):
    """ obtains the keywords defined in the input by the user
    """
    
    # Obtain the appropriate species string
    lvl_str = get_lvl_str(inp_str, lvl)

    # Get a list of the defined keywords
    defined_keys = []
    for line in lvl_str.splitlines():
        if '=' in line:
            tmp = line.strip().split('=')[0]
            defined_keys.append(tmp.strip())

    return defined_keys

def get_attr_call(key):
    calls = {'program' :  'get_prog', 
             'method'  :  'get_meth', 
             'basis'   :  'get_basis', 
              'ncycles':  'get_key_int',
              'mem'    :  'get_key_int',
              'nprocs' :  'get_key_int',
              'econv'  :  'get_key_str',
              'gconv'  :  'get_key_str',
              'orb_res'  :  'get_key_str'
       }
    if key in calls:
        return calls[key]
    else:
        log.warning('No method call for trait {}'.format(key))
        return 'error_call'

def error_call(s, lvl):
    log.warning('for theory {}'.format(lvl))
    return None

def get_key_int(inp_str, lvl, key):
    # Find the value for a given key, make int
    pattern = (key + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                      capturing(one_or_more(NONSPACE)))
    # Obtain the appropriate species string
    patt_str = get_lvl_str(inp_str, lvl)
    # Obtain the charge string
    patt_str = first_capture(pattern, patt_str)
    assert patt_str is not None
    return int(patt_str )

def get_key_str(inp_str, lvl, key):
    # Find the value for a given key, make int
    pattern = (key + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                      capturing(one_or_more(NONSPACE)))

    # Obtain the appropriate species string
    patt_str = get_lvl_str(inp_str, lvl)

    # Obtain the charge string
    patt_str = first_capture(pattern, patt_str)
    assert patt_str is not None
    return patt_str 

def get_prog(inp_str, lvl):
    # Find the value for a given key, make int
    pattern = ('program' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                      capturing(one_or_more(NONSPACE)))

    # Obtain the appropriate species string
    patt_str = get_lvl_str(inp_str, lvl)

    # Obtain the charge string
    patt_str = first_capture(pattern, patt_str)
    assert patt_str is not None
    if 'gaussian' in patt_str:
        patt_str = 'g09'
    return patt_str 

def get_meth(inp_str, lvl):
    # Find the value for a given key, make int
    pattern = ('method' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                      capturing(one_or_more(NONSPACE)))

    # Obtain the appropriate species string
    patt_str = get_lvl_str(inp_str, lvl)

    # Obtain the charge string
    patt_str = first_capture(pattern, patt_str)
    assert patt_str is not None
    return patt_str 
 
def get_basis(inp_str, lvl):
    # Find the value for a given key, make int
    pattern = ('basis' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                      capturing(one_or_more(NONSPACE)))

    # Obtain the appropriate species string
    patt_str = get_lvl_str(inp_str, lvl)

    # Obtain the charge string
    patt_str = first_capture(pattern, patt_str)
    assert patt_str is not None
    return patt_str

# OBTAIN PARAMETERS FOR ELSTRUCT JOBS

def get_ireact(species_file_name, level):
    """ Get the ireact parameter used to define transition states
    """

    # Set the ireact pattern
    ireact_pattern = ('program' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                     capturing(INTEGER))

    # Obtain the appropriate species string
    species_str = get_species_str(species_file_name, species)

    # Obtain the charge string
    ireact_str = first_capture(ireact_pattern, species_str)

    # Set the ireact parameter as an integer
    assert ireact_str is not None
    ireact = int(ireact_str)

    return ireact


