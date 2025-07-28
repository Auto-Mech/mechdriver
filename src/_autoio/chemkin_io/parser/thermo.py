""" functions operating on the thermo block string
"""

import itertools
import numpy as np
import autoparse.pattern as app
import autoparse.find as apf


COMMENTS_PATTERN = app.escape('!') + app.capturing(
    app.one_or_more(app.WILDCARD2))


def create_spc_nasa7_dct(block_str):
    """ Creates a spc_nasa7_dct


        :return spc_nasa7_dct: dictionary with spc names
                as keys and NASA-7 info as values

    """
    # Get the list of entries; at least four lines each
    # (possibly more for comment lines)
    entry_lst = create_entry_list(block_str)
    # Get the default midpoint temp
    default_temp_limits = get_default_temp_limits(block_str)
    default_midpoint = default_temp_limits[1]
    many_default_midpoints = list(
        itertools.repeat(default_midpoint, times=len(entry_lst)))
    # creates an iterator

    # Get the spc names, which will be the keys of the dictionary
    spc_names = list(map(get_spc_name, entry_lst))

    # Get the NASA-7 parameters, which will be the values of the dictionary
    nasa7_params = list(
        zip(
            map(get_notes, entry_lst),
            map(get_composition, entry_lst),
            map(get_phase, entry_lst),
            map(get_temp_limits, entry_lst, many_default_midpoints),
            map(get_coeffs, entry_lst)
        )
    )

    spc_nasa7_dct = dict(zip(spc_names, nasa7_params))

    return spc_nasa7_dct


# def create_entry_list(block_str, add_spaces=True):
def create_entry_list(block_str):
    """ Creates a list with each line of the thermo block_str as an
        entry in that list

        :param block_str: string for thermo block

        :return line_lst: list of strs for each line

    """
    block_str = apf.remove(COMMENTS_PATTERN, block_str)
    line_lst = list(apf.split_lines(block_str))

    # Get the indices of the entry header lines
    header_idxs = []
    for idx, line in enumerate(line_lst):
        if len(line) >= 80 and line[79] == '1':
            header_idxs.append(idx)
    if not header_idxs:
        raise ImportError(
            'No thermo headers in the file could be read.'
            + ' A "1" in column 80 marks this.')
    header_idxs.append(len(line_lst))
    # adds an entry to mark the end of the block_str

    # Check that there are at least 4 lines in each entry
    for idx, difference in enumerate(np.diff(header_idxs)):
        assert difference >= 4, (
            'There are less than 4 lines for this thermo entry: \n'
            f'{line_lst[header_idxs[idx]]:f}'
        )

    # Put together the thermo entries into a list of tuples
    entry_lst = []
    for idx in range(len(header_idxs)-1):  # skips the last entry
        entry = line_lst[header_idxs[idx]:header_idxs[idx+1]]
        entry_lst.append(entry)

    # Fix the problem of the missing spaces at
    # the beginning of lines without negative signs
    entry_lst = fix_lines(entry_lst)

    return entry_lst


def get_spc_name(entry, print_notice=False):
    """ Get the species name from the entry

        :param entry: the 4 or more lines making up the thermo entry
        :type entry: tuple
        :return spc_name: the name of the species, contained
                in the first 16 characters of the first line
        :rtype: str
    """
    first_line = entry[0]
    spc_name = first_line[0:25]  # first 25 characters (names should be 16)
    spc_name = spc_name.split()[0]  # just get the first space-separated entry
    if len(spc_name) > 16:
        print(f"Name '{spc_name}' is longer than 16 chars, "
              "this could cause parsing issues if a comment is included")

    # If there's something other than a species name here, print notice of this
    if len(spc_name) > 1 and print_notice:
        print(
            'In the first 16 characters of the following line, there is' +
            ' something other than only a species name.')
        print('(These extra characters will not be saved)')
        print(f'{first_line}')

    return spc_name


def get_notes(entry):
    """ Get the misc notes from the first line of the entry.
        Note that this is technically for the date,
        but it can contain whatever information is so desired.

        :param entry: the 4 or more lines making up the thermo entry
        :type entry: tuple
        :return notes: the freeform notes from the first line of the entry
        :rtype: str
    """
    first_line = entry[0]
    notes = first_line[18:24]  # characters 19 through 24

    return notes


def get_composition(entry):
    """ Get the molecular composition information from the first line of the entry


        :return composition: a dictionary containing the molecular composition
        :rtype: dct {element: quantity, ...}

    """
    first_line = entry[0]
    full_comp_str = first_line[24:44]  # characters 25 through 44

    # Will add the ability to read the actual values later
    # ...for now, just return the str
    # for idx in range(4):
    #     start = idx*5
    #     end = start + 5
    #     comp_str = full_comp_str[start:end]
    #     elem = comp_str[0:2]
    #     count = comp_str[2:5]

    composition = full_comp_str  # need to change later!

    # Also, there is room for a fifth entry right before the '1' on this line.
    # Need to add this later.

    return composition


def get_temp_limits(entry, default_midpoint):
    """ Get the temperatures from the first line of the entry


        :return temp_limits
        :rtype: list[floats] [low_limit, high_limit, midpoint]
                note the odd order
    """
    first_line = entry[0]
    formatted_entry = reform_entry(entry)  # for error printing
    # full_temp_str = first_line[45:73]  # characters 46 through 73

    # Read the high and low limits
    try:
        low_limit = float(first_line[45:55])  # characters 46 through 55
        high_limit = float(first_line[55:65])  # characters 56 through 65
    except ValueError:
        fline = first_line.split()
        low_limit = float(fline[-4])
        high_limit = float(fline[-3])
        midpoint = float(fline[-2])
        print(
            'Error processing the high and/or low temperatures' +
            f' in the following entry:\n{formatted_entry}')
        fline = first_line.split()
        low_limit = float(fline[-4])
        high_limit = float(fline[-3])
        midpoint = float(fline[-2])

    # If the midpoint read fails, replace it with the default value
    try:
        midpoint = float(first_line[65:73])  # characters 66 through 73
    except ValueError as valerr:
        if default_midpoint is None:
            raise ImportError(
                'No default mipoint temp in the file and no midpoint temp' +
                f' for the following entry:\n{formatted_entry}'
            ) from valerr
        midpoint = default_midpoint
        print(
            f'Using default midpoint temperature, {default_midpoint}' +
            f', for the following entry:\n{formatted_entry}'
        )

    temp_limits = [low_limit, high_limit, midpoint]

    return temp_limits


def get_phase(entry):
    """ Get the phase from the first line of the entry

        :return phase: 'G', 'L', or 'S'
        :rtype: str
    """
    first_line = entry[0]
    phase = first_line[44]  # character 45

    return phase


def get_coeffs(entry):
    """ Get the NASA-7 polynomial coefficients from
        the last three lines of the entry

        :return coeffs: the two sets of 7 polynomial coefficients
        :rtype: tuple(lsts)  ([high_coeffs], [low_coeffs])
    """
    formatted_entry = reform_entry(entry)
    line_counter = 1
    for idx, line in enumerate(entry):
        # If on the first line or if the line is too short, skip the line
        if idx == 0 or len(line) < 80:
            continue

        # If on the second line of the actual thermo data
        if line_counter == 1:
            try:
                high_coeffs = list((
                    float(line[0:15]), float(line[15:30]), float(line[30:45]),
                    float(line[45:60]), float(line[60:75])))
            except ValueError:
                print(
                    f'Error reading values in line {idx+1}' +
                    f' of the following entry:\n{formatted_entry}')
            line_counter += 1

        # If on the third line of the actual thermo data
        elif line_counter == 2:
            try:
                high_coeffs.extend(list((
                    float(line[0:15]), float(line[15:30]))))
                low_coeffs = list((
                    float(line[30:45]), float(line[45:60]),
                    float(line[60:75])))
            except ValueError:
                print(
                    f'Error reading values in line {idx + 1}' +
                    f' of the following entry:\n{formatted_entry}')
            line_counter += 1

        # If on the fourth line of the actual thermo data
        elif line_counter == 3:
            try:
                low_coeffs.extend(list((
                    float(line[0:15]), float(line[15:30]), float(line[30:45]),
                    float(line[45:60]))))
            except ValueError:
                print(
                    f'Error reading values in line {idx + 1}' +
                    f' of the following entry:\n{formatted_entry}')

    # Make sure three lines were read
    assert line_counter == 3, (
        'Less than three lines of coefficients were read' +
        f' for the following entry:\n{formatted_entry}'
    )
    coeffs = tuple((high_coeffs, low_coeffs))

    return coeffs


def get_default_temp_limits(block_str):
    """ Gets the default temperatures from the header of a thermo block str

    """
    block_str = apf.remove(COMMENTS_PATTERN, block_str)
    line_lst = list(apf.split_lines(block_str))

    # Loop over each line
    for line in line_lst:
        try:
            # Note that this order is different than
            # that in the individual thermo entries
            low_limit = float(line[0:10])
            midpoint = float(line[10:20])
            high_limit = float(line[20:30])
            break
        except ValueError:
            pass
        # Check if the first thermo entry has been reached
        if len(line) >= 80 and line[79] == '1':
            low_limit, midpoint, high_limit = None, None, None

    # Note that this order is different than
    # that in the individual thermo entries
    default_temp_limits = [low_limit, midpoint, high_limit]

    return default_temp_limits


def reform_entry(entry):
    """ Put the entry back together for error printing purposes

    """
    for idx, line in enumerate(entry):
        if idx == 0:
            formatted_entry = line
        else:
            formatted_entry += '\n' + line

    return formatted_entry


def fix_lines(entry_lst):
    """ Because Python is dumb and gets rid of leading whitespaces after a newline
    """
    for idx1, entry in enumerate(entry_lst):
        for idx2, line in enumerate(entry):
            if line:  # if the line is not empty
                first_char = line[0]
                # check if the first character is a digit
                # or a dot
                if first_char.isdigit() or first_char == '.':
                    line = ' ' + line
                if len(line) < 80: # for lines that start with more than one blank character
                    line = ' '*(80 - len(line)) + line
                entry[idx2] = line
        entry_lst[idx1] = entry

    return entry_lst
