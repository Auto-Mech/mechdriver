"""
library of reader functions for the species file
"""

import sys
from ast import literal_eval
from autoparse.find import first_capture
from autoparse.find import all_captures
from autoparse.pattern import capturing
from autoparse.pattern import zero_or_more
from autoparse.pattern import one_or_more
from autoparse.pattern import one_of_these
from autoparse.pattern import NONNEWLINE
from autoparse.pattern import NONSPACE
from autoparse.pattern import SPACE
from autoparse.pattern import WILDCARD
from autoparse.pattern import INTEGER
from autoparse.pattern import FLOAT
from qcelemental import constants as qcc
   
import logging
log   = logging.getLogger(__name__)


# Set the required keywords
REQUIRED_KEYWORDS = [
    'mult', 
    'charge', 
]

# Set the supported keywords for minima
SUPPORTED_KEYWORDS = [  
  'charge',
  'mult',     
  'geom',   
  'mc_nsamp',          
  'mc_tau',    
  'hind_inc',    
  'elec_levels',   
  'sym_factor',
  'inchi',
  'smiles',
  'geom'
]

# Set the supported keywords for transition states
SUPPORTED_KEYWORDS_TS = SUPPORTED_KEYWORDS + [
'isite',
'jsite',   
'ksite',   
'ireact',  
'aabs1',  
'aabs2',  
'babs1',  
'babs2',  
'babs3'  
]

CALLS = {    'geom'  :    'get_geom', 
             'inchi' :    'get_ichii', 
             'smiles':    'get_smiles', 
             'charge':    'get_charge', 
             'mult'  :    'get_mult', 
             'mc_nsamp':  'get_mc_nsamp', 
             'mc_tau'  :  'get_mc_tau',
             'hind_inc':  'get_hind_inc',
             'elec_levels':'get_elec_levels',
             'sym_factor': 'get_symmetry_factor'
       }
# FUNCTION TO READ IN A STRING FOR A SPECIFIC SPECIES # 

def get_species_str(inp_str, species):
    """ species input
    """
    # Set the species string pattern
    species_pattern = ('species' + one_or_more(SPACE) + species +  
                      capturing(one_or_more(WILDCARD, greedy=False)) +
                      'end')

    # Obtain the species string
    species_str = first_capture(species_pattern, inp_str)

    # make whole string lowercase?

    # test print
    #print(species_str)
    #sys.exit()

    return species_str


# FUNCTIONS TO CHECK WHICH SPECIES AND KEYWORDS HAVE BEEN DEFINED

def get_defined_species(inp_str):
    """ gets a list which specifies what species have been defined
    """
    species_def_pattern =  ('species' + one_or_more(SPACE) + capturing(one_or_more(NONSPACE)))  

    defined_species = all_captures(species_def_pattern, inp_str)

    return defined_species


def get_defined_keywords(inp_str, species):
    """ obtains the keywords defined in the input by the user
    """
    
    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Get a list of the defined keywords
    defined_keys = []
    for line in species_str.splitlines():
        if '=' in line:
            tmp = line.strip().split('=')[0]
            defined_keys.append(tmp.strip())

    return defined_keys


def check_required_keywords(inp_str, species):
    """ checks if certain keywords required by the program are defined in the species section 
    """
    
    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)
    
    # Obtain the defined keywords
    defined_keywords = get_defined_keywords(inp_str, species)

    # Verify all requird keywords are present
    required_isdefined = set(REQUIRED_KEYWORDS).issubset(set(defined_keywords))

    # Get the list of required keywords that are not currently defined 

    return required_isdefined 


def check_supported_keywords(inp_str, species):
    """ checks if keywords are defined that are not supported by the program
    """
    
    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)
    
    # Obtain the defined keywords
    defined_keywords = get_defined_keywords(inp_str, species)

    # Verify all required keywords are present
    keys_issupported = set(defined_keywords).issubset(set(SUPPORTED_KEYWORDS))
    
    # Get the keywords defined that are not currently supported

    return keys_issupported 


# FUNCTIONS TO READ PARAMETERS FOR MINIMA AND TRANSITION STATES #
def get_attr_call(key, repress = False):
    calls = CALLS
    if key in calls:
        return calls[key]
    else: 
        if not repress:
            log.warning('No method call for trait {}'.format(key))
        return 'error_call'

def error_call(s, spc):
    log.warning('for species {}'.format(spc))
    return None
    
def get_geom(inp_str, species):
    """ Get the electric charge of the species
    """

    # Set the geom pattern
    geom_pattern = ('geom' + zero_or_more(SPACE) + zero_or_more(SPACE) + "="  + 
                    capturing(one_or_more(WILDCARD, greedy=False)) +
                    one_of_these(list(CALLS.keys())))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the geom string
    geom_str = first_capture(geom_pattern, species_str)

    # Set the geom as an integer
    #assert geom_str is not None
    geom = geom_str
    return geom

def get_inchi(inp_str, species):
    """ Get the inchi name of the species
    """

    # Set the geom pattern
    ich_pattern = ('inchi' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                    capturing(one_or_more(NONNEWLINE)))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    #assert geom_str is not None
    return first_capture(ich_pattern, species_str)

def get_smiles(inp_str, species):
    """ Get the inchi name of the species
    """

    # Set the geom pattern
    smi_pattern = ('smiles' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                    capturing(one_or_more(NONNEWLINE)))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    #assert geom_str is not None
    return first_capture(smi_pattern, species_str)

def get_hind_inc(inp_str, species):
    """ Get the electric charge of the species
    """

    # Set the charge pattern
    hind_tau_pattern = ('hind_tau' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                      capturing(one_or_more(NONNEWLINE)))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the geom string
    hind_tau_str = first_capture(hind_tau_pattern, species_str)

    # Set the geom as an integer
    if hind_tau_str is not None:
        hind_tau =   float(hind_tau_str)
        hind_tau =  hind_tau * qcc.conversion_factor('degree', 'radian')

    else:
        hind_tau = 360. * qcc.conversion_factor('degree', 'radian')

    return hind_tau


def get_mc_tau(inp_str, species):
    """ Get the electric charge of the species
    """

    # Set the charge pattern
    mc_tau_pattern = ('mc_tau' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                      '"' +
                      capturing(one_or_more(WILDCARD, greedy=False)) +
                      '"')

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the geom string
    mc_tau_str = first_capture(mc_tau_pattern, species_str)

    # Set the geom as an integer
    if mc_tau_str is not None:
        if mc_tau_str == 'auto':
            mc_tau_dict = "auto"
        else:
            mc_tau_dict = literal_eval(mc_tau_str)
    else:
        mc_tau_dict = {}

    return mc_tau_dict


def get_charge(inp_str, species):
    """ Get the spin chargeiplicity of the species
    """

    # Set the charge pattern
    charge_pattern = ('charge' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                      capturing(INTEGER))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the charge string
    charge_str = first_capture(charge_pattern, species_str)

    # Set the charge as an integer
    assert charge_str is not None
    charge = int(charge_str)

    return charge


def get_mult(inp_str, species):
    """ Get the spin multiplicity of the species
    """

    # Set the mult pattern
    mult_pattern = ('mult' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                    capturing(INTEGER))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the charge string
    mult_str = first_capture(mult_pattern, species_str)

    # Set the charge as an integer
    assert mult_str is not None
    mult = int(mult_str)

    return mult


def get_mc_nsamp(inp_str, species):
    """ Get the number of samples for the Monte Carlo sampling
    """

    # Set the nsamp pattern
    mc_nsamp_pattern = ('mc_nsamp' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                        capturing(one_or_more(NONNEWLINE)))
    
    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the charge string
    mc_nsamp_str = first_capture(mc_nsamp_pattern, species_str)

    # Set the nsamp (maybe check the ntau variable)
    if mc_nsamp_str is not None:
        mc_nsamp_str = mc_nsamp_str.split() 
        if len(mc_nsamp_str) < 4:
            mc_nsamp = [False, 0, 0, 0, 0, int(mc_nsamp_str[0])]
        else:
            mc_nsamp = [True]
            for n in mc_nsamp_str:
                mc_nsamp.append(int(n))
            mc_nsamp.append(12)
    else:
        mc_nsamp = 0
    return mc_nsamp


def get_symmetry_factor(inp_str, species):
    """ read the number of samples for the Monte Carlo sampling
    """

    # Set the sym_factor pattern
    sym_factor_pattern = ('sym_factor' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                        capturing(FLOAT))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the charge string
    sym_factor_str = first_capture(sym_factor_pattern, species_str)

    # Set the sym_factor as an float
    if sym_factor_str is not None:
        sym_factor = float(sym_factor_str)
    else:
        sym_factor = 1.00

    return sym_factor


def get_elec_levels(inp_str, species):
    """ Get the electronic levels of the species
    """

    # Set the electronic energy levels pattern (as a list-of-lists)
    elec_levels_pattern = ('elec_levels' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                      '"' +
                      capturing(one_or_more(WILDCARD, greedy=False)) +
                      '"')

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the geom string
    elec_levels_str = first_capture(elec_levels_pattern, species_str)

    # Set the geom as an integer
    if elec_levels_str is not None:
        elec_levels_list = literal_eval(elec_levels_str)
    else:
        elec_levels_list = [[1.0, 0.0]]

    return elec_levels_list


# FUNCTIONS TO READ PARAMETERS FOR THE TRANSITION STATE # 

def get_isite(inp_str, species):
    """ Get the isite parameter used to define transition states
    """

    # Set the isite pattern
    isite_pattern = ('isite' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                     capturing(INTEGER))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the charge string
    isite_str = first_capture(isite_pattern, species_str)

    # Set the isite parameter as an integer
    assert isite_str is not None
    isite = int(isite_str)

    return isite


def get_jsite(inp_str, species):
    """ Get the jsite parameter used to define transition states
    """

    # Set the jsite pattern
    jsite_pattern = ('jsite' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                     capturing(INTEGER))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the charge string
    jsite_str = first_capture(jsite_pattern, species_str)

    # Set the jsite parameter as an integer
    assert jsite_str is not None
    jsite = int(jsite_str)

    return jsite


def get_ksite(inp_str, species):
    """ Get the ksite parameter used to define transition states
    """

    # Set the ksite pattern
    ksite_pattern = ('ksite' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                     capturing(INTEGER))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the charge string
    ksite_str = first_capture(ksite_pattern, species_str)

    # Set the ksite parameter as an integer
    assert ksite_str is not None
    ksite = int(ksite_str)

    return ksite


def get_ireact(inp_str, species):
    """ Get the ireact parameter used to define transition states
    """

    # Set the ireact pattern
    ireact_pattern = ('ireact' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                     capturing(INTEGER))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the charge string
    ireact_str = first_capture(ireact_pattern, species_str)

    # Set the ireact parameter as an integer
    assert ireact_str is not None
    ireact = int(ireact_str)

    return ireact


def get_aabs1(inp_str, species):
    """ Get the aabs1 parameter used to define transition states
    """

    # Set the aabs1 pattern
    aabs1_pattern = ('aabs1' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                     capturing(FLOAT))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the charge string
    aabs1_str = first_capture(aabs1_pattern, species_str)

    # Set the aabs1 parameter as an integer
    assert aabs1_str is not None
    aabs1 = float(aabs1_str)

    return aabs1


def get_aabs2(inp_str, species):
    """ Get the aabs2 parameter used to define transition states
    """

    # Set the aabs2 pattern
    aabs2_pattern = ('aabs2' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                     capturing(FLOAT))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the charge string
    aabs2_str = first_capture(aabs2_pattern, species_str)

    # Set the aabs2 parameter as an integer
    assert aabs2_str is not None
    aabs2 = float(aabs2_str)

    return aabs2


def get_babs1(inp_str, species):
    """ Get the babs1 parameter used to define transition states
    """

    # Set the babs1 pattern
    babs1_pattern = ('babs1' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                     capturing(FLOAT))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the charge string
    babs1_str = first_capture(babs1_pattern, species_str)

    # Set the babs1 parameter as an integer
    assert babs1_str is not None
    babs1 = float(babs1_str)

    return babs1


def get_babs2(inp_str, species):
    """ Get the babs2 parameter used to define transition states
    """

    # Set the babs2 pattern
    babs2_pattern = ('babs2' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                     capturing(FLOAT))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the charge string
    babs2_str = first_capture(babs2_pattern, species_str)

    # Set the babs2 parameter as an integer
    assert babs2_str is not None
    babs2 = float(babs2_str)

    return babs2


def get_babs3(inp_str, species):
    """ Get the babs3 parameter used to define transition states
    """

    # Set the babs3 pattern
    babs3_pattern = ('babs3' + zero_or_more(SPACE) + '=' + zero_or_more(SPACE) + 
                     capturing(FLOAT))

    # Obtain the appropriate species string
    species_str = get_species_str(inp_str, species)

    # Obtain the charge string
    babs3_str = first_capture(babs3_pattern, species_str)

    # Set the babs3 parameter as an integer
    assert babs3_str is not None
    babs3 = float(babs3_str)

    return babs3
