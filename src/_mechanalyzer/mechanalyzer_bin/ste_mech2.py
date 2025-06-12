""" Script to generate files with stereoexpanded
    reaction and species list for a mechanism.
"""

import os
import argparse
import time
import itertools as it

import ioformat
import automol
import chemkin_io
import mechanalyzer
from mechanalyzer.builder import sorter


def main(workdir,
         run_expansion=True,
         run_reduction='enant',
         check_mechanism=False,
         nprocs='auto'):
    """ carry out all the mechanism things you wanna do
    """

    # Check
    print(f'Running Expansion: {run_expansion}')
    print(f'Running Reduction: {run_reduction}')
    print('---------\n\n')

    # Read and parse the build and mechanism file
    print('\n---- Parsing the species and mechanism files ---\n')

    bld_str = ioformat.pathtools.read_file(
        workdir, 'build.dat',
        remove_comments='#', remove_whitespace=True)
    file_dct, _ = mechanalyzer.parser.build_input_file(bld_str)
    (rxn_param_dct, mech_spc_dct,
     isolate_spc, sort_lst) = input_from_working_dictionary(workdir, file_dct)

    # Run algorithms that were requested
    sccs_rxn_dct_lst, ccs_sccs_spc_dct = None, None

    if run_expansion:
        # Expansion mechanism, generating mech,spc files for all S-CCSs
        sccs_rxn_dct_lst, ccs_sccs_spc_dct = expand_stereo(
            rxn_param_dct, mech_spc_dct,
            check_mechanism=check_mechanism,
            nprocs=nprocs)

        # Write the files of all the generated S-CCSs
        write_all_sccs_mechfiles(
            sccs_rxn_dct_lst,
            workdir, file_dct['out_spc'], file_dct['out_mech'],
            sort_lst, isolate_spc)

        for idx, lst in enumerate(sccs_rxn_dct_lst):
            print(idx)
            print(lst)
            print('----')

    if run_reduction:

        # Read in all of the mechanism files of the expaneded S-CCSs, ifneeded
        if not run_expansion:
            print('Reading expanded mech files for spc and reaction lists')
            ccs_sccs_spc_dct, name_ich_dct = read_all_sccs_spc(
                workdir, file_dct['out_spc'])
            sccs_rxn_dct_lst = read_all_sccs_rxns(
                workdir, file_dct['out_mech'], name_ich_dct)

        for idx, lst in enumerate(sccs_rxn_dct_lst):
            print(idx)
            print(lst)
            print('----')

        # Select the (CCS, S-CCS) surfaces for the reduced mechanism
        best_combo, best_combo_ichs = reduction_selection(
            ccs_sccs_spc_dct, algorithm=run_reduction)

        # Add in any required diastereomer surfaces precluded by the reduction
        print('combo test 1', best_combo)
        best_combo += find_diastereomer_abstraction_sccs(
            sccs_rxn_dct_lst, ccs_sccs_spc_dct,
            best_combo, best_combo_ichs)
        print('combo test 2', best_combo)

        # Write the final, reduced mechanism
        write_best_combination(
            best_combo, workdir,
            file_dct['out_mech'], file_dct['out_spc'], file_dct['sort'])


# EXPANSION ALGORITHMS #
def expand_stereo(rxn_param_dct, mech_spc_dct,
                  check_mechanism=False, nprocs='auto'):
    """ Generate the fully expanded mechanism
    """

    # Remove reactions that should not be there
    if check_mechanism:
        print('\n Removing improper reactions')
        rxn_param_dct = mechanalyzer.builder.remove_improper_reactions(
            rxn_param_dct, mech_spc_dct)

    # Expand the stereo of the species and reactions
    print('\n---- Adding stereochemistry to InChIs of mechanism'
          ' species where needed ---\n')
    mech_spc_dct = mechanalyzer.parser.spc.stereochemical_spc_dct(
        mech_spc_dct, nprocs=nprocs, all_stereo=True)
    print('Mechanism species with stereo added')
    for name, dct in mech_spc_dct.items():
        print(f'Name: {name:<25s} InChI: {dct["inchi"]}')

    print('\n---- Expanding the list of mechanism reactions to include all'
          ' valid, stereoselective permutations ---\n')
    sccs_rxn_dct_lst = mechanalyzer.builder.expand_mech_stereo(
        rxn_param_dct, mech_spc_dct, nprocs=nprocs)
    ccs_sccs_spc_dct = _build_ccs_sccs_spc_dct(sccs_rxn_dct_lst)

    return sccs_rxn_dct_lst, ccs_sccs_spc_dct


# REDUCTION ALGORITHMS #
def reduction_selection(ccs_sccs_spc_dct, algorithm='enant'):
    """ Determine the (CCS, S-CCS) pairs that best represent a sufficient
        mechanism with all required stereoisomers and minimizes redundant
        enantiomers.
    """

    if algorithm == 'enant':
        sccs_combo_lst = _all_sccs_combos(ccs_sccs_spc_dct)
        best_combo, best_combo_ichs = _reduce_via_enantiomers(
            ccs_sccs_spc_dct, sccs_combo_lst)
    elif algorithm == 'overlap':
        best_combo, best_combo_ichs = _reduce_via_max_overlap(
            ccs_sccs_spc_dct)

    return best_combo, best_combo_ichs


def _reduce_via_enantiomers(ccs_sccs_spc_dct, combo_lst):
    """ find the sccs combo that has the shortest unique spc lst
        which means it has the least number of enantiomers
        reduce by doing all possible combinations of sccss and figuring out
        which combination results in the shortest unique species list
    """
    uniq_spc_lst_lst = ()
    num_uniq_spc_for_combo = ()
    print('Checking {:g} combinations'.format(len(combo_lst)))
    for combo in combo_lst:
        spc_ich_lst = ()
        for idxs in combo:
            spc_ich_lst += ccs_sccs_spc_dct[idxs[0]][idxs[1]]
        uniq_spc_lst = ()
        for spc in spc_ich_lst:
            if spc not in uniq_spc_lst:
                uniq_spc_lst += (spc,)
        uniq_spc_lst_lst += (uniq_spc_lst,)
        num_uniq_spc_for_combo += (len(uniq_spc_lst),)
        # enant_count = 0
        # for spc_a, spc_b in it.combinations(spc_ich_lst, 2):
        #     if automol.chi.are_enantiomers(spc_a, spc_b):
        #         enant_count += 1
        # combo_ent_count += (enant_count,)
        # print('found {:g} enantiomers for this combo'.format(enant_count))
    # print(combo_ent_count)
    # min_enants = min(combo_ent_count)
    # print('minimum overlap', min_enants)
    # best_combo_idx = combo_ent_count.index(min_enants)
    # best_combo = combo_ent_count[best_combo_idx]
    # print('best combo', best_combo)
    shortest_spc_lst = min(num_uniq_spc_for_combo)
    print('greatest overlap', shortest_spc_lst)
    best_combo_idx = num_uniq_spc_for_combo.index(shortest_spc_lst)
    best_combo = combo_lst[best_combo_idx]
    print('best combo', best_combo)
    enant_count = 0
    for spc_a, spc_b in it.combinations(uniq_spc_lst_lst[best_combo_idx], 2):
        if automol.chi.are_enantiomers(spc_a, spc_b):
            print('Enantiomer pair:', spc_a, spc_b)
            enant_count += 1
    print('found {:g} enantiomers for this combo'.format(enant_count))

    best_combo_ichs = ()
    for ccs, sccs in best_combo:
        best_combo_ichs += ccs_sccs_spc_dct[ccs][sccs]

    return best_combo, best_combo_ichs


def _group_sccs(ccs_sccs_spc_dct):
    """ using the keys of the ccs_sccs_spc_dct group all
        ccs combos into a tuple of ccs_sccs tuples
    """
    all_idxs = ()
    ccs_lst = sorted(list(ccs_sccs_spc_dct.keys()))
    for ccs in ccs_lst:
        sccs_lst = sorted(list(ccs_sccs_spc_dct[ccs].keys()))
        sccs_tup = tuple((ccs, sccs) for sccs in sccs_lst)
        all_idxs += (sccs_tup,)
    return all_idxs


def _separate_single_ccs(all_idxs):
    sccs_lsts = ()
    ccs_lsts = ()
    for idx_lst in all_idxs:
        if len(idx_lst) == 1:
            ccs_lsts += idx_lst
        else:
            sccs_lsts += (idx_lst,)
    return sccs_lsts, ccs_lsts


def _all_sccs_combos(ccs_sccs_spc_dct):
    """ build all possible combinations of sccs
    """
    combo_lst = ((),)
    all_idxs = _group_sccs(ccs_sccs_spc_dct)
    sccs_lsts, ccs_lsts = _separate_single_ccs(all_idxs)
    for ccs_tups in sccs_lsts:
        new_combo_lst = ()
        print('Building full combo list from', ccs_tups)
        for idxs in ccs_tups:
            for combo in combo_lst:
                new_combo_lst += (combo + (idxs,),)
        combo_lst = new_combo_lst

    new_combo_lst = ()
    for combo in combo_lst:
        new_combo_lst += (combo + ccs_lsts,)
    return new_combo_lst


def _reduce_via_max_overlap(ccs_sccs_spc_dct):
    """ Find the best combination of CCS,S-CCS to make the reduced
        mechanism by choosing an S-CCS that overlaps best with the previosly
        chosen mechanism. For the first CCS, we choose the first S-CCS
        arbitrarily as a starting point
    """

    # Choose the sccs for each ccs that works with the rest of the mech
    best_combo = ()
    best_combo_ichs = ()
    for ccs_idx, sccs_dct in ccs_sccs_spc_dct.items():
        if not best_combo:
            best_combo += ((ccs_idx, 0),)
            best_combo_ichs += tuple(sccs_dct[0])
            continue
        max_overlap = 0
        best_idx = 0
        for sccs_idx, spc_dct in sccs_dct.items():
            num_overlap = len(set(best_combo_ichs) & set(spc_dct))
            if num_overlap > max_overlap:
                max_overlap = num_overlap
                best_idx = sccs_idx
        best_combo += ((ccs_idx, best_idx),)
        best_combo_ichs += tuple(sccs_dct[best_idx])

    return best_combo, best_combo_ichs


# FUNCTIONS TO REPAIR REDUCED LISTS #
def find_diastereomer_abstraction_sccs(sccs_rxn_dct_lst, ccs_sccs_spc_dct,
                                       best_combo, best_combo_ichs):
    """ Grab any missed diatereomer S-CCSs
    """

    print('Identifying diastereomer abstractions')
    dias_idxs = mechanalyzer.builder.diastereomer_abstractions(
        sccs_rxn_dct_lst, ccs_sccs_spc_dct,
        best_combo,
        best_combo_ichs)

    if dias_idxs:
        dia_str = ', '.join(
            (f'({ccs}, {sccs})' for (ccs, sccs) in dias_idxs)
        )
        print(f'Found new diastereomers to add: {dia_str}')

    return dias_idxs


# HELPERS #
def _dictionaries_from_rxn_lst(sccs_rxn_lst):
    """ transform the reaction list to the dictionaries the writer likes
    """
    ste_mech_spc_dct, ste_rxn_dct = {}, {}
    ste_mech_spc_dct = mechanalyzer.builder.update_spc_dct_from_reactions(
        sccs_rxn_lst, ste_mech_spc_dct)
    ste_rxn_dct = mechanalyzer.builder.update_rxn_dct(
        sccs_rxn_lst, ste_rxn_dct, ste_mech_spc_dct)
    ste_mech_spc_dct = mechanalyzer.builder.remove_spc_not_in_reactions(
        ste_rxn_dct, ste_mech_spc_dct)
    return ste_mech_spc_dct, ste_rxn_dct


def _build_ccs_sccs_spc_dct(sccs_rxn_dct_lst):
    """ Build the species dictionary for each (CCS, S-CCS) surface
        from its reactions
    """

    ccs_sccs_spc_dct = {}
    for ccs_idx, sccs_rxn_dct in enumerate(sccs_rxn_dct_lst):
        for sccs_idx, sccs_rxn_lst in sccs_rxn_dct.items():
            ste_mech_spc_dct, _ = _dictionaries_from_rxn_lst(
                sccs_rxn_lst)
            is_valid = mechanalyzer.builder.valid_enantiomerically(
                ste_mech_spc_dct)
            if is_valid:
                spc_ichs = tuple(dct['inchi']
                                 for dct in ste_mech_spc_dct.values())
                if ccs_idx not in ccs_sccs_spc_dct:
                    ccs_sccs_spc_dct[ccs_idx] = {0: spc_ichs}
                else:
                    ccs_sccs_spc_dct[ccs_idx][sccs_idx] = spc_ichs

    return ccs_sccs_spc_dct


# I/O #
def input_from_working_dictionary(workdir, file_dct):
    """ read mechanism info using file dictionary
        and current working directory
    """

    inp_spc_str, inp_mech_str, sort_str = _input_info_from_file(
        workdir, file_dct['inp_spc'], file_dct['inp_mech'], file_dct['sort'])

    mech_spc_dct = mechanalyzer.parser.spc.build_spc_dct(inp_spc_str, 'csv')
    rxn_param_dct = mechanalyzer.parser.mech.parse_mechanism(
        inp_mech_str, 'chemkin')
    isolate_spc, sort_lst = mechanalyzer.parser.mech.parse_sort(sort_str)

    return rxn_param_dct, mech_spc_dct, isolate_spc, sort_lst


def _input_info_from_file(workdir, species_name, mech_name, sort_name):
    """ read mechanism info using filenames
    """
    inp_spc_str = ioformat.pathtools.read_file(
        workdir, species_name,
        remove_comments='!', remove_whitespace=True)
    inp_mech_str = ioformat.pathtools.read_file(
        workdir, mech_name,
        remove_comments='#', remove_whitespace=True)
    sort_str = ioformat.pathtools.read_file(
        workdir, sort_name,
        remove_comments='#', remove_whitespace=True)

    return inp_spc_str, inp_mech_str, sort_str


def read_all_sccs_spc(workdir, spc_file_name):
    """ create a nested dictionary with ccs, sccs, spc_lst
        where spc lst is the inchis of that ccs, sccs combo
    """
    files = os.listdir(workdir)
    spc_files = [fil for fil in files if spc_file_name + '_' in fil]

    ccs_sccs_spc_dct = {}
    name_ich_dct = {}
    for fil in spc_files:
        ccs, sccs = [int(idx) for idx in fil.split('_')[-2:]]
        spc_lines = io.read_file(fil).splitlines()
        name_col = [
            idx for idx, val in enumerate(spc_lines[0].split("'"))
            if val == 'name'][0]
        ich_col = [
            idx for idx, val in enumerate(spc_lines[0].split("'"))
            if val == 'inchi'][0]

        spc_names = tuple(line.split("'")[name_col] for line in spc_lines[1:])
        spc_ichs = tuple(line.split("'")[ich_col] for line in spc_lines[1:])
        if ccs not in ccs_sccs_spc_dct:
            ccs_sccs_spc_dct[ccs] = {sccs: spc_ichs}
        else:
            ccs_sccs_spc_dct[ccs][sccs] = spc_ichs

        for name, ich in zip(spc_names, spc_ichs):
            name_ich_dct.update({name: ich})

    return ccs_sccs_spc_dct, name_ich_dct


def read_all_sccs_rxns(workdir, mech_file_name, name_ich_dct):
    """ created nested dicinary with ccs, sccs where the list of reactions
        are the values
    """
    files = os.listdir(workdir)
    mech_files = [fil for fil in files if mech_file_name + '_' in fil]

    idxs = []
    for fil in mech_files:
        ccs, sccs = [int(idx) for idx in fil.split('_')[-2:]]
        idxs.append([ccs, sccs])

    sccs_rxn_dct = {}
    for ccs, sccs in sorted(idxs):
        fil = f'{mech_file_name}_{ccs}_{sccs}'
        _mech_str = ioformat.pathtools.read_file(
            workdir, fil,
            remove_comments='!', remove_whitespace=True)
        rxn_param_dct = mechanalyzer.parser.mech.parse_mechanism(
            _mech_str, 'chemkin')

        rxn_lst = ()
        for rxn in rxn_param_dct:
            rxn_lst += ((tuple(name_ich_dct[rct] for rct in rxn[0]),
                         tuple(name_ich_dct[prd] for prd in rxn[1])),)

        if ccs in sccs_rxn_dct:
            sccs_rxn_dct[ccs][sccs] = rxn_lst
        else:
            sccs_rxn_dct[ccs] = {sccs: rxn_lst}

    sccs_rxn_dct_lst = []
    for ccs, sccs_dct in sccs_rxn_dct.items():
        _dct = {}
        for sccs, rxn_lst in sccs_dct.items():
            _dct.update({sccs: rxn_lst})
            print('read mechfile', ccs, sccs)
        sccs_rxn_dct_lst.append(_dct)

    return sccs_rxn_dct_lst


def write_all_sccs_mechfiles(sccs_rxn_dct_lst,
                             workdir, out_spc_name, out_mech_name,
                             sort_lst, isolate_spc):
    """ write the files for each S-CCS
    """

    full_ste_mech_spc_dct, full_ste_rxn_dct = {}, {}
    for ccs_idx, sccs_rxn_dct in enumerate(sccs_rxn_dct_lst):
        for sccs_idx, sccs_rxn_lst in sccs_rxn_dct.items():

            # Get rxns and spc of each individual S-CCS
            ste_mech_spc_dct, ste_rxn_dct = _dictionaries_from_rxn_lst(
                sccs_rxn_lst)
            is_valid = mechanalyzer.builder.valid_enantiomerically(
                ste_mech_spc_dct)

            # Write mechanism files of S-CCS
            if is_valid:
                _write_mechanism(
                    ste_mech_spc_dct, ste_rxn_dct,
                    workdir, out_spc_name, out_mech_name,
                    sort_lst.copy(), isolate_spc, ccs_idx, sccs_idx)

                # Combine all rxn and spc into full dictionary
                full_ste_mech_spc_dct.update(ste_mech_spc_dct)
                full_ste_rxn_dct.update(ste_rxn_dct)

    # Write mechanism files for fully expanded mechanism with all stereoisomers
    _write_mechanism(
        full_ste_mech_spc_dct, full_ste_rxn_dct,
        workdir, out_spc_name, out_mech_name,
        sort_lst.copy(), isolate_spc, 9999, 9999)


def write_best_combination(best_combo, workdir,
                           out_mech_name, out_spc_name, sort_name):
    """ make combined species and mechanism dictionaries
        from a set of sccs
    """

    mech_spc_dct = {}
    rxn_param_dct = {}
    for (ccs, sccs) in best_combo:
        print('writing from', ccs, sccs)
        mech_file = '{}_{:g}_{:g}'.format(out_mech_name, ccs, sccs)
        spc_file = '{}_{:g}_{:g}'.format(out_spc_name, ccs, sccs)
        inp_info = _input_info_from_file(
            workdir, spc_file, mech_file, sort_name)
        next_spc_dct = mechanalyzer.parser.spc.build_spc_dct(
            inp_info[0], 'csv')
        next_rxn_param_dct = mechanalyzer.parser.mech.parse_mechanism(
            inp_info[1], 'chemkin')
        mech_spc_dct.update(next_spc_dct)
        rxn_param_dct.update(next_rxn_param_dct)
        for rxn in next_rxn_param_dct:
            print(f'  {rxn}')

    isolate_spc, sort_lst = mechanalyzer.parser.mech.parse_sort(inp_info[2])

    _write_mechanism(
        mech_spc_dct, rxn_param_dct,
        workdir, 'red_species.csv', 'red_mechanism.dat',
        sort_lst, isolate_spc)


def _write_mechanism(ste_mech_spc_dct, ste_rxn_dct,
                     workdir, out_spc_name, out_mech_name,
                     sort_lst, isolate_spc,
                     ccs_idx=None, sccs_idx=None):
    """ write mechanism
    """

    # Write the new (sorted) species dictionary to a string
    headers = ('smiles', 'inchi', 'inchikey', 'mult', 'charge')
    ste_mech_spc_dct_sort = mechanalyzer.parser.spc.reorder_by_atomcount(
        ste_mech_spc_dct)
    csv_str = mechanalyzer.parser.spc.csv_string(
        ste_mech_spc_dct_sort, headers)

    # Write initial string to call the sorter
    mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
        elem_tuple=None,
        mech_spc_dct=ste_mech_spc_dct_sort,
        spc_nasa7_dct=None,
        rxn_param_dct=ste_rxn_dct,
        rxn_cmts_dct=None)

    # Use strings to generate ordered objects
    param_dct_sort, ste_mech_spc_dct_sort, cmts_dct, _, _ = sorter.sorted_mech(
        csv_str, mech_str, isolate_spc, sort_lst)
    rxn_cmts_dct = chemkin_io.writer.comments.get_rxn_cmts_dct(
        rxn_sort_dct=cmts_dct)

    # Write the dictionaries to ordered strings
    csv_str = mechanalyzer.parser.spc.csv_string(
        ste_mech_spc_dct_sort, headers)
    mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
        elem_tuple=(),
        mech_spc_dct=ste_mech_spc_dct_sort,
        spc_nasa7_dct=None,
        rxn_param_dct=param_dct_sort,
        rxn_cmts_dct=rxn_cmts_dct)

    # Write the species and mechanism files
    if ccs_idx is not None and sccs_idx is not None:
        ioformat.pathtools.write_file(
            csv_str, workdir,
            out_spc_name + '_{:g}_{:g}'.format(ccs_idx, sccs_idx))
        ioformat.pathtools.write_file(
            mech_str, workdir,
            out_mech_name + '_{:g}_{:g}'.format(ccs_idx, sccs_idx))
    else:
        ioformat.pathtools.write_file(csv_str, workdir, out_spc_name)
        ioformat.pathtools.write_file(mech_str, workdir, out_mech_name)


if __name__ == '__main__':

    # Check input
    PAR = argparse.ArgumentParser()
    PAR.add_argument(
        '-e', '--expansion', default=True, type=bool,
        help='Generate stereoexpanded surfaces [True(def), False]')
    PAR.add_argument(
        '-r', '--reduction', default='enant', type=str,
        help='Use reduction algorithm [enant(def), overlap]')
    PAR.add_argument(
        '-n', '--nprocs', default='auto', type=str,
        help='Number of processors to use [auto(def)]')
    OPTS = vars(PAR.parse_args())

    # Convert number of processors to an int if an integer is given
    if OPTS['nprocs'] != 'auto':
        try:
            OPTS['nprocs'] = int(OPTS['nprocs'])
        except ValueError:
            print('nprocs not set to auto or an integer')
            sys.exit()

    # Initialize the start time for script execution
    t0 = time.time()

    # Set useful global variables
    CWD = os.getcwd()

    # Execute all desired algorithms
    main(CWD,
         run_expansion=OPTS['expansion'],
         run_reduction=OPTS['reduction'],
         check_mechanism=False,
         nprocs=OPTS['nprocs'])

    # Compute script run time and print to screen
    tf = time.time()
    print('\n\nScript executed successfully.')
    print(f'Time to complete: {tf-t0:.2f}')
    print('Exiting...')
