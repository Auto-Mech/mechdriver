"""
Parses a spc.csv file containing species information to obtain a mech_spc_dct
"""

import csv
import copy
from ioformat import pathtools
import automol
from automol.chi import inchi_to_amchi
from automol.chi import smiles as ich_to_smi
from automol.chi import canonical_enantiomer
from automol.chi import formula as ich_to_fml
from automol.chi import low_spin_multiplicity as _low_spin_mult
from automol.form import from_string as str_to_fml

ALLOWED_COL_NAMES = (
    'name',
    'smiles',
    'inchi',
    'inchikey',
    'mult',
    'charge',
    'exc_flag',
    'fml',
    'sens',
    'hof',
    'canon_enant_ich',
    'fct_grp'
)

TRIP_DCT = {
    'O':        'InChI=1S/O',
    'O2':       'InChI=1S/O2/c1-2',
    'CH2':      'InChI=1S/CH2/h1H2',
    'C3H2':     'InChI=1S/C3H2/c1-3-2/h1-2H',
    'H2CCC':    'InChI=1S/C3H2/c1-3-2/h1H2',
    'C6H6O-B1': 'InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-2,4-5H,3H2',
    'CH2CHN':   'InChI=1S/C2H3N/c1-2-3/h2H,1H2',
}

CMTS = '!'  # character used to define comments


def load_mech_spc_dcts(filenames, path, quotechar="'", chk_ste=False,
                       chk_match=False, verbose=True, canon_ent=True, allow_multiple=False):
    """ Obtains multiple mech_spc_dcts given a list of spc.csv filenames

        :param filenames: filenames of the spc.csv files to be read
        :type filename: list
        :param quotechar: character used to optionally ignore commas; " or '
        :type quotechar: str
        :param chk_ste: whether or not to check inchis for stereo completeness
        :type chk_ste: Bool
        :param chk_match: whether or not to check that inchis and smiles match
        :type chk_match: Bool
        :param verbose: whether or not to print lots of warnings
        :type verbose: Bool
        :return mech_spc_dcts: list of mech_spc_dcts
        :rtype: list
    """

    mech_spc_dcts = []
    for filename in filenames:
        print(f'Loading mech_spc_dct for the file {filename}...')
        mech_spc_dct = load_mech_spc_dct(
            filename, path, quotechar=quotechar, chk_ste=chk_ste,
            chk_match=chk_match, verbose=verbose, canon_ent=canon_ent,
            allow_multiple=allow_multiple)
        mech_spc_dcts.append(mech_spc_dct)

    return mech_spc_dcts


def load_mech_spc_dct(filename, path, quotechar="'", chk_ste=False,
                      chk_match=False, verbose=True, canon_ent=True, allow_multiple=False):
    """ Obtains a single mech_spc_dct given a spc.csv filename

        :param filename: filename of the spc.csv file to be read
        :type filename: str
        :param quotechar: character used to optionally ignore commas; " or '
        :type quotechar: str
        :param chk_ste: whether or not to check inchis for stereo completeness
        :type chk_ste: Bool
        :param chk_match: whether or not to check that inchis and smiles match
        :type chk_match: Bool
        :param verbose: whether or not to print lots of warnings
        :type verbose: Bool
        :param canon_ent: add/check if canonical enantiomer
        :type canon ent: Bool
        :param allow_multiple: add multiple occurrences in spc dct as lumped entries
        :type allow_multiple: Bool
        :return mech_spc_dct: identifying information on species in a mech
        :rtype: dct {spc1: spc_dct1, spc2: ...}
    """

    file_str = pathtools.read_file(path, filename, print_debug=True)
    mech_spc_dct = parse_mech_spc_dct(
        file_str, quotechar=quotechar, chk_ste=chk_ste, chk_match=chk_match,
        verbose=verbose, canon_ent=canon_ent, allow_multiple=allow_multiple)

    return mech_spc_dct


def parse_mech_spc_dct(file_str, quotechar="'", chk_ste=False,
                       chk_match=False, verbose=True, canon_ent=True, allow_multiple=False):
    """ Obtains a single mech_spc_dct given a string parsed from a spc.csv file

        :param file_str: the string that was read directly from the .csv file
        :type file_str: str
        :param quotechar: the quotechar used to optionally ignore commas
        :type quotechar: str
        :param chk_ste: whether or not to check inchis for stereo completeness
        :type chk_ste: Bool
        :param chk_match: whether or not to check that inchis and smiles match
        :type chk_match: Bool
        :param canon_ent: add/check if canonical enantiomer
        :type canon ent: Bool
        :param allow_multiple:  add multiple occurrences in spc dct as lumped entries
        :type allow_multiple: Bool
        :param verbose: whether or not to print lots of warnings
        :type verbose: Bool
        :return mech_spc_dct: identifying information on species in a mech
        :rtype: dct {spc1: spc_dct1, spc2: ...}
    """

    # Remove comment lines
    file_str = pathtools.remove_comment_lines(file_str, CMTS)

    # Check for incorrect quotechar usage
    if quotechar == '"':
        assert "'" not in file_str, (
            'A quotechar input of double quotation marks was selected, but at '
            'least one instance of a single quotation mark exists in the '
            'csv file. Use only double quotation marks.')
    elif quotechar == "'":
        assert '"' not in file_str, (
            'A quotechar input of single quotation marks was selected, but at '
            'least one instance of a double quotation mark exists in the '
            'csv file. Use only single quotation marks.')
    else:
        raise NotImplementedError(
            f'The quotechar {quotechar} is not allowed. Options are either '
            f'a single or double quotation mark.')

    # Split into lines
    lines = file_str.split('\n')
    if lines[-1] == '':
        del lines[-1]  # gets rid of the annoying last line that gets added

    # Build the mech_spc_dct line by line
    mech_spc_dct = {}
    dup_spc_idx = {} # indexes for duplicate species
    found_inchis = []
    for idx, line in enumerate(lines):
        if idx == 0:
            headers = parse_first_line(line, quotechar=quotechar)
            if canon_ent and 'canon_enant_ich' not in headers:
                print("Warning: user selected the 'canon_ent' option, but the"
                      " field 'canon_enant_ich' is not in the CSV file.\n"
                      "The canonical enantiomer will have to be calculated "
                      "for every species, which might be slow.")
        else:
            cols = parse_line(line, idx, headers, quotechar=quotechar)
            if cols is not None:
                spc, spc_dct = make_spc_dct(cols, headers)
                # Fill in the spc_dct and then add it to the mech_spc_dct
                spc_dct = fill_spc_dct(
                    spc_dct, spc, chk_ste=chk_ste,
                    chk_match=chk_match, canon_ent=canon_ent)
                # Check that the species name was not already defined
                if allow_multiple:
                    if spc in mech_spc_dct:
                        # skip if InChI is the same- probably just written twice by mistake
                        if spc_dct['inchi'] in found_inchis:
                            continue

                        if spc in dup_spc_idx:
                            lumped_idx = dup_spc_idx[spc] + 1
                        else:
                            lumped_idx = 1

                        dup_spc_idx[spc] = lumped_idx
                        lumped_spc_name = '{}-LUMPED-{}'.format(spc, lumped_idx)

                        # otherwise assign as lumped
                        print('The species {} was already identified. '
                            'The first and current InChIs are \n'
                            '{} \n'
                            '{} \n'
                              'Assuming lumped species. Adding current InChI as {}'.format(
                                spc, spc_dct['inchi'], mech_spc_dct[spc]['inchi'],
                                lumped_spc_name)
                        )
                        spc = lumped_spc_name
                else:
                    assert spc not in mech_spc_dct, (
                        f'The species name {spc} appears in the csv file more than'
                        f' once. The second time is on line {idx + 1}, {line}.')

                mech_spc_dct[spc] = spc_dct
                found_inchis.append(spc_dct['inchi'])
            else:
                print(f'Line {idx + 1} appears to be empty. Skipping...')

    # Find species with the same chemical identifiers but different names
    check_for_dups(mech_spc_dct, printwarnings=verbose)

    return mech_spc_dct


def make_spc_dct(cols, headers):
    """ Makes a spc_dct by reading parsed contents of a single line

        :param cols: a list of the entries in the line
        :type cols: list
        :param headers: column headers (i.e., first line of the csv)
        :type headers: list
        :return spc: the mechanism name for the species
        :rtype: str
        :return spc_dct: identifying information for a single species
        :rtype: dct
    """

    # Build the spc_dct for this species, column by column
    spc_dct = {}
    for idx, col in enumerate(cols):
        header = headers[idx]
        if header == 'name':
            spc = col
        # If charge or mult or exc_flag, convert to int if not empty
        elif header in ('charge', 'mult', 'exc_flag') and col != '':
            spc_dct[header] = int(col)
        # If on formula, convert to dict
        elif header == 'fml':
            spc_dct[header] = str_to_fml(col)
            # if 'fct_grp' is an header: do not try to add
        elif header == 'fct_grp':
            # if not dictionary: convert
            spc_dct[header] = fct_grp_todct(
                    col)
        else:
            spc_dct[header] = col

    # If inchi or smiles not in headers, fill in blank value to avoid KeyError
    # (can't both be missing; see parse_first_line)
    if 'inchi' not in headers:
        spc_dct['inchi'] = ''
    elif 'smiles' not in headers:
        spc_dct['smiles'] = ''

    return spc, spc_dct


def fill_spc_dct(spc_dct, spc, chk_ste=True, chk_match=True, canon_ent=True):
    """ Fills in missing values in a spc_dct

        :param spc_dct: identifying information for a single species
        :type: dct
        :param spc: species name
        :type spc: str
        :param chk_ste: whether or not to check inchis for stereo completeness
        :type chk_ste: Bool
        :param chk_match: whether or not to check that inchis and smiles match
        :type chk_match: Bool
        :return full_spc_dct: beefed-up spc_dct
        :rtype: dct
    """
    # print('Filling species dictionary with missing values (InChI/SMILES, AMChI, stereo if selected). \
    #     Consider preprocessing with mechanalyzer_bin/expand_species.py first if not done')

    full_spc_dct = copy.deepcopy(spc_dct)

    # Add SMILES or InChI if missing
    assert full_spc_dct['inchi'] != '' or full_spc_dct['smiles'] != '', (
        f'For {spc}, both inchi and smiles are empty')

    if full_spc_dct['inchi'] == '':
        smi = full_spc_dct['smiles']
        check_smi(smi, spc)
        full_spc_dct['inchi'] = automol.smiles.chi(smi)

    elif full_spc_dct['smiles'] == '':
        ich = full_spc_dct['inchi']
        check_ich(ich, spc, chk_ste=chk_ste)
        full_spc_dct['smiles'] = ich_to_smi(ich)

    else:  # if smiles and inchi were both already included
        if chk_match:
            check_smi_and_ich(
                full_spc_dct['smiles'], full_spc_dct['inchi'], spc,
                chk_ste=chk_ste)
        else:
            check_ich(full_spc_dct['inchi'], spc, chk_ste=chk_ste)

    # add AMChI
    full_spc_dct = mech_inchi_to_amchi({spc: full_spc_dct}, convert=canon_ent)[spc]

    if not canon_ent:
        full_spc_dct['canon_enant_ich'] = full_spc_dct['inchi']
    elif 'canon_enant_ich' not in full_spc_dct:
        full_spc_dct = add_canonical_enantiomer(
            {spc: full_spc_dct})[spc]

    # Add charge and exc_flag if missing; assume 0 for both
    if 'charge' not in full_spc_dct or full_spc_dct['charge'] == '':
        full_spc_dct['charge'] = 0
    if 'exc_flag' not in full_spc_dct or full_spc_dct['exc_flag'] == '':
        full_spc_dct['exc_flag'] = 0

    # Add multiplicity if missing by guessing
    if 'mult' not in full_spc_dct or full_spc_dct['mult'] == '':
        ich = full_spc_dct['inchi']
        smi = full_spc_dct['smiles']
        mult = _low_spin_mult(ich)
        if mult == 1 and ich in TRIP_DCT.values():
            sing = False
            if '(S)' in smi or '(S)' in ich or '(S)' in spc:
                sing = True
            elif 'singlet' in smi or 'singlet' in ich or 'singlet' in spc:
                sing = True
            if not sing:  # otherwise, just leave as 1
                mult = 3
        full_spc_dct['mult'] = mult

    # Add formula if missing
    if 'fml' not in full_spc_dct or full_spc_dct['fml'] == '':
        fml = ich_to_fml(full_spc_dct['inchi'])
        full_spc_dct['fml'] = fml

    return full_spc_dct


def parse_first_line(first_line, quotechar="'"):
    """ Parses the first line in the spc.csv file

        :param first_line: a string with the contents of the first line in the
            csv file
        :type first_line: str
        :param quotechar: the quotechar used to optionally ignore commas
        :type quotechar: str
        :return headers: column headers (i.e., first line of the csv)
        :rtype: list
    """

    # Remove the UTF encoding if present
    if '\ufeff' in first_line:
        first_line = first_line.replace('\ufeff', '')

    # Parse the line
    headers = next(csv.reader([first_line], quotechar=quotechar))

    # Check for appropriateness of column headers (case insensitive)
    for header in headers:
        assert header.lower() in ALLOWED_COL_NAMES, (
            f'Header {header} is not allowed. Options: {ALLOWED_COL_NAMES}.')
    assert 'name' in headers, (
        "The species name, 'name', must be included in the csv file headers.")
    assert 'inchi' in headers or 'smiles' in headers, (
        "At least one of the following chemical identifiers must be included"
        f" in the csv file headers: 'inchi' or 'smiles'. Headers: {headers}")

    return headers


def parse_line(line, idx, headers, quotechar="'"):
    """ Parses a line in the spc.csv file (other than the first line)

        :param line: a string with the contents of a single line in the csv
            file
        :type line: str
        :param idx: the index of the line
        :type idx: int
        :param headers: column headers (i.e., first line of the csv)
        :type headers: list
        :param quotechar: the quotechar used to optionally ignore commas
        :type quotechar: str
        :return cols: a list of the entries in the line
        :rtype: list
    """
    # The `csv` module can't handle trailing spaces, so strip the line before reading
    cols = next(csv.reader([line.strip()], quotechar=quotechar))
    if line == '':
        cols = None
    else:
        assert len(cols) == len(headers), (
            f'Line {idx + 1}, {line}, has {len(cols)} columns. The first line '
            f'in the csv file indicates {len(headers)} columns are needed.')

    return cols


def check_ich(ich, spc, chk_ste=True):
    """ Checks an inchi for various problems

        :param ich: inchi string
        :type ich: str
        :param spc: species name
        :type spc: str
        :param chk_ste: whether or not to check inchi for stereo completeness
        :type chk_ste: Bool
        :return error: whether or not an error was found with the inchi
        :rtype: Bool
    """

    error = False
    gra = automol.chi.graph(ich)
    if gra is None:
        print(f"The spc '{spc}' has an invalid InChI, '{ich}'")
        error = True

    # If indicated, check for stereochemical completeness of inchi
    if chk_ste:
        if not automol.chi.is_complete(ich):
            print(f"The spc '{spc}' has a stereo-invalid InChI, '{ich}'")
            error = True

    assert not error, (
        'Errors when checking InChI validity')


def check_smi(smi, spc):
    """ Checks a smiles for various problems

        :param smi: smiles string
        :type smi: str
        :param spc: species name
        :type spc: str
        :return error: whether or not an error was found with the smiles
        :rtype: Bool
    """

    error = False
    gra = automol.smiles.graph(smi)
    if gra is None:
        print(f"The spc '{spc}' has an invalid SMILES string, '{smi}'")
        error = True

    assert not error, (
        'Errors in SMILES dictionary!')



def check_smi_and_ich(smi, ich, spc, chk_ste=True):
    """ Checks a smiles and inchi for various problems

        :param smi: smiles string
        :type smi: str
        :param ich: inchi string
        :type ich: str
        :param spc: species name
        :type spc: str
        :param chk_ste: whether or not to check inchi for stereo completeness
        :type chk_ste: Bool
        :return error: whether or not an error was found with inchi or smiles
        :rtype: Bool
    """

    # Check the smiles and inchi for correctness
    error1 = check_smi(smi, spc)
    error2 = check_ich(ich, spc, chk_ste=chk_ste)

    # Recalculate inchi from smiles and see if it's the same
    # Note: only do if both initial tests passed and if not checking stereo
    error3 = False
    if not any([error1, error2]) and not chk_ste:
        recalc_ich = automol.smiles.chi(smi)
        if ich != recalc_ich:
            print(f"The spc '{spc}' has non-matching SMILES & InChI strings")
            error3 = True

    error = any([error1, error2, error3])

    assert not error, (
        'Errors in checking smiles VS InChIs')


def check_for_dups(mech_spc_dct, printwarnings=True):
    """ Checks a mech_spc_dct for species that are identical except in name
    """

    def are_spcs_same(spc_dct1, spc_dct2):
        """ Checks if two species are chemically identical
        """

        ich1 = spc_dct1['inchi']
        mlt1 = spc_dct1['mult']
        chg1 = spc_dct1['charge']
        exc1 = spc_dct1['exc_flag']

        ich2 = spc_dct2['inchi']
        mlt2 = spc_dct2['mult']
        chg2 = spc_dct2['charge']
        exc2 = spc_dct2['exc_flag']

        same = all([ich1 == ich2, mlt1 == mlt2, chg1 == chg2, exc1 == exc2])

        return same

    spcs = tuple(mech_spc_dct.keys())  # convert to tuple for indexing
    spc_dcts = tuple(mech_spc_dct.values())
    for outer_idx, outer_spc in enumerate(spcs):
        outer_dct = spc_dcts[outer_idx]
        for inner_idx in range(outer_idx + 1, len(spcs)):
            inner_spc = spcs[inner_idx]
            inner_dct = spc_dcts[inner_idx]
            if are_spcs_same(outer_dct, inner_dct) and printwarnings:
                print(f'{outer_spc} and {inner_spc} are chemical twins!')

def mech_inchi_to_amchi(mech_spc_dct, convert = True):
    """ convert inchi to amchi where needed """
    for spc_dct in mech_spc_dct.values():
        try:
            spc_dct['inchi'] = inchi_to_amchi(spc_dct['inchi'], convert=convert)
        except: #if errors are raised
            print('*Warning- failed conversion to AMChI, but spc will be read')
            continue

    return mech_spc_dct

def add_canonical_enantiomer(mech_spc_dct, dummy=False):
    """ add canonical enantiomer to species dictionaries
    """
    for spc_dct in mech_spc_dct.values():
        if 'canon_enant_ich' not in spc_dct:
            if not dummy:
                spc_dct['canon_enant_ich'] = canonical_enantiomer(
                    spc_dct['inchi'])
            else:
                spc_dct['canon_enant_ich'] = spc_dct['inchi']
    return mech_spc_dct

def add_fct_grp_dct(mech_spc_dct, species_subset = 'all'):
    """ add functional group dictionary to species dictionary
        mech_spc_dct needs to have either smiles or inchi to get the graph
    """
    for spc, dct in mech_spc_dct.items():

        # if there is a specific subset of species requested: filter
        if isinstance(species_subset, list):
            if spc not in species_subset:
                mech_spc_dct[spc]['fct_grp'] = {}
                continue
        # now check geometry
        if 'inchi' in dct.keys():
            ich = dct['inchi']
        elif 'smiles' in dct.keys():
            ich = automol.smiles.chi(dct['smiles'])
        try:
            gra = automol.geom.graph(automol.chi.geometry(ich))
        #except AssertionError as e:
        except:
            #print(f'{spc} skipped for failure in geom generation {e}')
            print(f'{spc} skipped for failure in geom generation')
            mech_spc_dct[spc]['fct_grp'] = {}
            continue

        mech_spc_dct[spc]['fct_grp'] = automol.graph.classify_species(gra)

    return mech_spc_dct


def fct_grp_tostr(fct_grp_dct):
    """ turns functional group dictionary into a string
        {'A1-M':1, 'C5-RSR': 2} => 'A1-M:1 C5-RSR:2'
    """
    return ' '.join(f"{k}:{v}" for k, v in fct_grp_dct.items())

def fct_grp_todct(fct_grp_str):
    """ turns functional group string into a dictionary
        'A1-M:1 C5-RSR:2' => {'A1-M':1, 'C5-RSR': 2}
    """
    return {k: int(v) for k, v in (item.split(':') for item in fct_grp_str.split())}
