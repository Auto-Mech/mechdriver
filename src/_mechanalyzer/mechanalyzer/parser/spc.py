""" Handles the construction, parsing, and manipulation of
    mechanism species dictionaries structured as

    {name: spc_data_dct}

    Also handles ancillary dictionaries used contain and
    organize information related to the species dictionary
"""

import os
import copy
import csv
import pandas as pd
import automol
from automol.form import string
from autorun import timeout, execute_function_in_parallel
import ioformat.pathtools as text_parser
import thermfit
from mechanalyzer.parser.csv_ import csv_dct
from mechanalyzer.parser.new_spc import fct_grp_tostr

# LIST SETTING THE STANDARD ORDER OF HEADERS
STD_HEADERS = (
    'name',
    'smiles',
    'inchi',
    'inchikey',
    'mult',
    'charge'
)


# Write a spc_dct to a CSV string
def csv_string(spc_dct, headers):
    """ Convert a mechanism species dictionary to a CSV file formatted
        string. The headers provided correspond to the physical data
        or parameters that will be parsed from the dictionary and written
        into the file.

        Note that the first column of the CSV string will always
        be 'name' and so that should be precluded from the list. The
        remaining headers will be written in the order provided.

        :param spc_dct: mechanism species dictionary
        :type spc_dct: dict[str: dict]
        :param headers: headers of the CSV file
        :type headers: tuple(str)
        :rtype: str
    """

    # Build spc sub dct of just info in the headers
    _csv_dct = {}
    for name, dct in spc_dct.items():
        _csv_dct[name] = {}
        for header in headers:
            #  if dictionaries are found: turn them into strings (example: fml, fct_grp_dct)
            val = dct.get(header, None)
            if isinstance(val, dict):
                if header == 'fml':
                    val = string(val)
                elif header == 'fct_grp':
                    val = fct_grp_tostr(val)
            _csv_dct[name][header] = val
            

    # Build the datagrame and resultant CSV string
    dframe = pd.DataFrame.from_dict(_csv_dct, orient='index')
    _csv_str = dframe.to_csv(
        index_label='name',
        quotechar="'",
        quoting=csv.QUOTE_NONNUMERIC)

    return _csv_str


# headers function i will probably move inside the function above
def csv_headers(spc_dct, include_missing=False):
    """ Determine what the headers should be for writing a csv string dct
    """

    # Build spc sub dct of just info in the headers
    headers = []
    for dct in spc_dct.values():
        headers += list(dct.keys())

    # Turn to set to remove duplicates
    # headers = ('name',) + tuple(set(headers))
    headers = tuple(set(headers))

    # Sort the headers by the standard list
    headers = automol.util.sort_by_list(
        headers, STD_HEADERS, include_missing=include_missing)

    return headers


# Build spc_dct from file/string i/o
def load_spc_dcts(spc_csv_filenames, direc):
    """ Reads one or more spc.csv files (in the AutoMech format). Outputs a
        list of spc_dcts.

        :param spc_csv_filenames: filenames containing spc info
        :type spc_csv_filenames: list [filename1, filename2, ...]
        :param direc: directory with file(s) (all must be in same directory)
        :type direc: str
        :return spc_dcts: list of spc_dcts
        :rtype: list of dcts [spc_dct1, spc_dct2, ...]
    """
    spc_dcts = []
    for spc_csv_filename in spc_csv_filenames:
        print('spc csv test:', spc_csv_filename)
        spc_csv_str = text_parser.read_file(direc, spc_csv_filename)
        spc_dct = build_spc_dct(spc_csv_str, 'csv')
        spc_dcts.append(spc_dct)

    return spc_dcts


def build_spc_dct(spc_str, spc_type):
    """ Get a dictionary of all the input species
        indexed by InChi string
    """

    if spc_type == 'csv':
        spc_dct = csv_dct(spc_str)
    else:
        raise NotImplementedError

    return spc_dct


# Modify an existing spc_dct
def reorder_by_atomcount(spc_dct):
    """ Returns a species dictionary ordered by increasing N of atoms
        in the species. Useful for nicer mechanism writing
    """

    natom_df = pd.Series(index=list(spc_dct.keys()))
    for key in spc_dct.keys():
        ich = spc_dct[key]['inchi']
        fml_dct = automol.chi.formula(ich)
        natoms = automol.form.atom_count(fml_dct)
        natom_df[key] = natoms
    natom_df = natom_df.sort_values(ascending=True)

    # Rewrite spc_dct
    sorted_idx = list(natom_df.index)
    sorted_val = list(map(spc_dct.get, sorted_idx))
    spc_dct_ordered = dict(zip(sorted_idx, sorted_val))

    return spc_dct_ordered


def add_heat_of_formation_basis(spc_dct,
                                ref_schemes=('cbh0', 'cbh1', 'cbh2'),
                                nprocs='auto'):
    """ Adds all species required to form the basis set for requested CBHn
        heat-of-formation calculations for all species in the mechanism
        species dictionary.

        Function will loop over all species in the input dictionary,
        determine the basis and check if each basis currently in dictionary.

        If species is missing it is assigned a name:
        cbh{N}_{smiles} where smiles is the smiles string, and N is the
        lowest order cbh scheme for which a basis species was found.
    """

    def _filter_redundant_basis(cbh_ref_dct, ref_scheme, cbh_smiles):
        """ Get only spc of cbh_ref_dct that have smiles not already added
            to spc_dctfrom earlier CBHn schemes
        """
        new_dct = {}
        for dct in cbh_ref_dct.values():
            tempn_smiles = automol.chi.smiles(dct['inchi'])
            if tempn_smiles not in cbh_smiles:
                # Add to new dictionry
                tempn = ref_scheme + '_' + tempn_smiles
                new_dct[tempn] = dct
                new_dct[tempn]['smiles'] = tempn_smiles
                # Add smiles to new smiles bookkeeping list
                cbh_smiles.append(tempn_smiles)

        return new_dct, cbh_smiles

    # Find all references
    cbh_smiles = []
    for ref_scheme in ref_schemes:
        basis_dct = thermfit.prepare_basis(
            ref_scheme, spc_dct, tuple(spc_dct.keys()),
            nprocs=nprocs, print_log=True)
        uniref_dct = thermfit.unique_basis_species(basis_dct, spc_dct)
        nonred_dct, cbh_smiles = _filter_redundant_basis(
            uniref_dct, ref_scheme, cbh_smiles)
        spc_dct.update(nonred_dct)

    return spc_dct


def add_instability_products(mech_spc_dct, nprocs='auto', stereo=True):
    """ Add species related to the instability products
    """

    print('Finding unstable molecules:')
    all_instab_ichs = ()
    for name, dct in mech_spc_dct.items():
        ich = dct['inchi']
        instab_ichs = automol.reac.instability_product_inchis(
            ich, stereo=stereo)
        if instab_ichs is not None:
            print(f'Found instability for {name} = {ich}')
            print(f'- {instab_ichs}')
            all_instab_ichs += instab_ichs

    if all_instab_ichs:

        print('Adding unstablie species to species dictionary')

        all_instab_ichs = tuple(set(all_instab_ichs))

        _name_ich_dct = name_inchi_dct(mech_spc_dct)
        _ich_name_dct = {ich: name for name, ich in _name_ich_dct.items()}
        for ich in all_instab_ichs:
            _name = _ich_name_dct.get(ich)
            if _name is None:
                _name = f'instab_{automol.chi.smiles(ich)}'
                mech_spc_dct[_name] = thermfit.create_spec(ich)
                print(f'{_name} = {ich} being added to species dictonary')
            else:
                print(f'{_name} = {ich} already in species dictonary')

    return mech_spc_dct


def stereochemical_spc_dct(
        spc_dct, nprocs='auto', all_stereo=False, enant=True):
    """ read the species file in a .csv format and write a new one
        that has stero information
    """

    # Generate names and procs for parallel stereochem add'n
    init_names = list(spc_dct.keys())

    # Add stereo using multiple processes
    args = (spc_dct, all_stereo, enant)
    ste_dcts = execute_function_in_parallel(
        _add_stereo_to_dct, init_names, args, nprocs=nprocs)

    # Build and reoroder the stereo dct (only works if all_stereo=False)
    ste_spc_dct = {}
    for dct in ste_dcts:
        ste_spc_dct.update(dct)
    if not all_stereo:
        ste_spc_dct_ord = {name: ste_spc_dct[name] for name in init_names
                           if name in ste_spc_dct}
    else:
        ste_spc_dct_ord = copy.deepcopy(ste_spc_dct)

    return ste_spc_dct_ord


def _add_stereo_to_dct(init_dct, all_stereo, enant, names, output_queue):
    """ Builds a modified species dictionary for a set of names where
        each sub species dictionary contains an InChI string with
        stereochemical layers being added.

        The function can either add all of the stereochemical species or
        just pick one.
    """

    # Functions for calling automol for ring counts and stereo add'n
    # @timeout(30)
    # def _nrings(name, dct):
    #     """ Count the number of rings
    #     """
    #     try:
    #         nrings = len(automol.graph.rings(
    #             automol.chi.graph(dct['inchi'])))
    #     except:
    #         print('Cannot produce graph for {} '.format(name))
    #         nrings = 2000
    #     return nrings

    @timeout(2000)
    def _generate_stereo(name, ich, all_stereo=False):
        """ stereo
        """
        ret_ichs, worked = [ich], True
        # print('inchi test:', name, ich)
        # print('complete inchi test:', automol.chi.is_complete(ich))
        # print('add_stereo  inchi test:', automol.chi.add_stereo(ich))
        # print('expand_stereo inchi test:', automol.chi.expand_stereo(ich))
        if all_stereo:
            try:
                if not automol.chi.is_complete(ich):
                    ret_ichs = (
                         automol.chi.expand_stereo(ich, enant=enant))
            except:  # noqa: E722
                print(f'{name} timed out in stereo generation')
                worked = False
        else: 
            try:
                if not automol.chi.is_complete(ich):
                    ret_ichs = (
                        [automol.chi.add_stereo(ich)])
            except:  # noqa: E722
                ich_attempt = automol.chi.expand_stereo(ich, enant=enant)
                if len(ich_attempt) > 0:
                    ret_ichs = [ich_attempt[0]]
                else:   
                    worked = False
        # Loop over strings and convert stereo inchi to amchi, if needed
        # may not work perfectly since you rely on the inchi code
        ret_ichs = [automol.chi.inchi_to_amchi(ich) for ich in ret_ichs]
        
        return ret_ichs, worked

    # Assess the species the code is able to add stereochemistry to
    names_str = ', '.join(names)
    print(f'Processor {os.getpid()} will check species: {names_str}')

    good_names = []
    for name in names:
        # if _nrings(name, init_dct[name]) < 10:
        good_names.append(name)
        # else:
        #    print('{} has >10 rings; stereo not implemented'.format(name))

    # Construct a new dictionary with stereochemical inchi strings
    new_dct = {}
    good_names_str = ', '.join(good_names)
    print(f'Processor {os.getpid()} will check species: {good_names_str}')
    for name in good_names:

        # Generate ichs with stereo and hashes
        ich = init_dct[name]['inchi']
        print('working on ich:', name, ich)
        ste_ichs, success = _generate_stereo(name, ich, all_stereo=all_stereo)
        if not success:
            print('failed on ich:', name, ich)
            continue
        print('succeeded to produce full ich for:', name, ich, ste_ichs)

        # Add name and inchi info to string
        spc_dct_keys = tuple(key for key in init_dct[name].keys()
                             if key not in ('inchi', 'inchikey'))
        for idx, ste_ich in enumerate(ste_ichs):
            # name gen wrong, should use formula builder
            sname = name+f'({str(idx+1)})' if idx != 0 else name
            new_dct[sname] = {
                'inchi': ste_ich,
                'inchikey': automol.chi.inchi_key(ste_ich)
            }
            for key in spc_dct_keys:
                new_dct[sname][key] = init_dct[name][key]

    output_queue.put((new_dct,))
    print(f'Processor {os.getpid()} finished')


# HELPER FUNCTIONS
# Handle Formula Count Objects
def name_inchi_dct(spc_dct):
    """Generates inchis dictionary from species dictionary
    """
    ich_dct = {}
    for key in spc_dct.keys():
        if 'ts' not in key and 'global' not in key:
            ich_dct[key] = spc_dct[key]['inchi']

    return ich_dct


# add functions like adding the hashkey
def add_hashkey(spc_dct):
    """ Add the inchi hashkey
    """

    for dct in spc_dct.values():
        ich = dct.get('inchi')
        ick = automol.chi.inchi_key(ich) if ich is not None else None
        dct.update({'inchikey': ick})

    return spc_dct
