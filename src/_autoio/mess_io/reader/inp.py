""" Intended to be a start for parsing MESS input files
"""

import copy
import autoparse.pattern as app
import autoparse.find as apf
from ioformat import headlined_sections
from ioformat import remove_comment_lines
import mess_io.reader._mess_keys as KEYS
from autoparse import cast as ap_cast
from phydat import phycon

UNITS = app.escape('[') + app.one_or_more(app.NONSPACE) + app.escape(']')
UNITS_CAP = app.escape('[') + app.capturing(app.one_or_more(app.NONSPACE)) +\
    app.escape(']')


def parse_model_section(section, keywd):
    """ Parses a single section of a model and parses
    """

    if keywd == 'EnergyRelaxationFlags':
        pass
    elif keywd == 'CollisionFrequency':
        pass
    elif keywd == 'EnergyRelaxation':
        pass
    elif keywd == 'Well':
        well = parse_well(section)
        quit()
    elif keywd == 'Barrier':
        pass
    elif keywd == 'Bimolecular':
        pass
    else:
        print(f'The keyword {keywd} has not been implemented')


def parse_well(well_section):

    keywds = KEYS.WELL
    sections, matches, _ = get_sections(well_section, keywds)
    for idx, section in enumerate(sections):
        match = matches[idx]
        if match == 'EnergyRelaxation':
            pass  # get energy_relaxation object
        elif match == 'Species':
            parse_species(section)
        elif match == 'Escape':
            pass  # get Escape info (probably not an object)
        elif match == 'CollisionFrequency':
            pass  # get collision object
        elif match == 'WellExtensionCap[kcal/mol]':
            pass  # get WellCap info (probably not an object)
        else:
            print(f'The keyword {match} has not been implemented')

        # Add the object to the well object

#    return well

def parse_species(species_section):

    keywds = KEYS.SPECIES
    sections, matches, _ = get_sections(species_section, keywds)
    #print('species_section:\n', species_section)
    for idx, section in enumerate(sections):
        match = matches[idx]
        if match == 'RRHO':
            parse_rrho(section)
        elif match == 'Variational':
            pass  # get Variational object
        else:
            print(f'The keyword {match} has not been implemented')

        # Add to species object

    #return species


def parse_rrho(rrho_section):


    keywds = KEYS.RRHO
    sections, matches, _ = get_sections(rrho_section, keywds)
    #print('rrho_section:\n', rrho_section)
    for idx, section in enumerate(sections):
        #print(section)
        #print('rrho match:\n', match)
        match = matches[idx]
        if 'ElectronicLevels' in match:
            get_elec_lvls(section)
        elif 'ZeroEnergy' in match:
            get_zero_egy(section)
        elif match == 'Variational':
            pass  # get Variational object
        else:
            #print(f'The keyword {match} has not been implemented')
            pass
        # Add to species object

    #return species


    


def get_elec_lvls(section):

    entry_pat = app.capturing(app.NUMBER)
    header_pat = 'ElectronicLevels' + UNITS_CAP + app.SPACES + entry_pat 

    # Find the header line (might not be the first line)
    for idx, line in enumerate(section):
        line = _rm_all_comments(line)
        captures = apf.all_captures(header_pat, line)
        if captures is not None:
            units, num_lvls = ap_cast(captures)[0]
            break  # exit after finding the header line

    # Look for the info from the next few lines
    elec_lvls = []
    for new_idx in range(num_lvls):
        line = section[idx + new_idx + 1]
        line = _rm_all_comments(line)
        vals = get_mult_vals(line)
        elec_lvls.append(vals)
        
    return elec_lvls, units


# Helper functions
def get_sections(lines, keywds):
    """ Takes a list of lines and breaks them into groups based on keywords.
        Any lines before the first keyword are returned as the separate object
        front_matter.

        :param lines: lines of text to be grouped into sections
        :type lines: list [str1, str2, ...]
        :param keywds: keywords to be used for breaking into sections
        :type keywds: list [str1, str2, ...]
        :return sections: list of sections; each section is a list of lines
        :rtype: list [section1, section2, ...]
        :return matches: the matching keywords for each section
        :rtype: list [keyword1, keyword2, ...]
        :return front_matter:
        :rtype: list [str1, str2, ...]
    """

    sections = []
    matches = []
    front_matter = None
    start_idx = None
    for idx, line in enumerate(lines):
        line = _rm_all_comments(line)
        for keywd in keywds:
            keywd_regex = get_keywd_regex(keywd)
            if apf.first_capture(keywd_regex, line) is not None:
                if idx == 0:  # if first match on first line, no front_matter
                    pass
                elif start_idx is None:  # if first match is not on first line
                    front_matter = lines[:idx]  # all lines up to current
                else:
                    sections.append(lines[start_idx:idx])
                matches.append(keywd)
                start_idx = idx
                break  # don't keep looking for keywords on this line

    # Check if no keywords were found
    if start_idx is None:
        print(f'Warning: none of the below keywords were found\n\n{keywds}')

    # Add the last section
    sections.append(lines[start_idx:len(lines)])

    return sections, matches, front_matter


def _conv_egy(val, units):
    """ Converts an input value to units of kcal

    """
    if units == '1/cm':
        val = val * phycon.WAVEN2KCAL
    elif units == 'kJ/mol':
        val = val * phycon.KJ2KCAL
    elif units == 'eV': 
        raise NotImplemented  # I don't see eV anywhere in phycon...

    return val


def get_mult_vals(line):

    pattern = app.capturing(app.NUMBER) + app.one_of_these([app.SPACES, 
                                                            app.LINE_END])
    captures = apf.all_captures(pattern, line)
    vals = list(ap_cast(captures))

    return vals


def _rm_all_comments(input_str):
    """ Removes all comments since MESS uses '!' or '#' for comments

    """
    filtered_str = copy.deepcopy(input_str)
    filtered_str = remove_comment_lines(
        filtered_str, delim_pattern=app.escape('!'))
    filtered_str = remove_comment_lines(
        filtered_str, delim_pattern=app.escape('#'))

    return filtered_str


def _keywd_units_val(keywd, line):
    """ Finds a keyword, the keyword's units, and a single value in a line

    """
    
    pattern = app.ZSPACES + keywd + UNITS_CAP + app.SPACES + app.capturing(
        NUMBER)
    line = _rm_all_comments(line)
    units, val = ap_cast(apf.all_captures(pattern, line))[0]
    
    return units, val 
          

#def _keywd_val(keywd, line):

    

def get_keywd_regex(keywd):

    # Add backslashes to escape some characters
    keywd = keywd.replace('/', '\/')
    keywd = keywd.replace('[', '\[')
    keywd = keywd.replace(']', '\]')

    regex = app.ZSPACES + keywd + app.one_of_these([app.SPACE, app.LINE_END])

    return regex
