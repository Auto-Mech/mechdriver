"""
Read the mechanism file
"""

import mechanalyzer
from lib.amech_io.parser import ptt


MECH_INP = 'inp/mechanism.dat'


def build_pes_dct(job_path, mech_type,
                  spc_dct, run_obj_dct, sort_rxns=True):
    """Build the PES dct
    """

    # Read the string
    mech_str = ptt.read_inp_str(job_path, MECH_INP, remove_comments='!')

    # Build the total PES dct
    mech_info = mechanalyzer.parser.pes.read_mechanism_file(
        mech_str, mech_type, spc_dct, sort_rxns=sort_rxns)
    pes_dct = mechanalyzer.parser.pes.build_pes_dct(*mech_info[1:])

    # Build an index dct relating idx to formula
    idx_dct, form_dct = build_pes_idx_dct(pes_dct)

    # Reduce the PES dct to only what the user requests
    if run_obj_dct:
        pesnums = [idx_pair[0] for idx_pair in run_obj_dct]
    else:
        pesnums = [idx_pair[0] for idx_pair in run_obj_dct]
    reduced_pes_dct = reduce_pes_dct_to_user_inp(pes_dct, pesnums)

    # Get a dct for all of the connected channels with the PESs to run
    conn_chnls_dct = mechanalyzer.parser.pes.connected_channels_dct(
        reduced_pes_dct)

    # Form the pes dct that has info formatted to run
    # Get the models in here
    run_pes_dct = pes_dct_w_rxn_lsts(
        reduced_pes_dct, idx_dct, form_dct, conn_chnls_dct, run_obj_dct)

    # Print the channels for the whole mechanism file
    print_pes_channels(pes_dct)

    return run_pes_dct


# FUNTIONS FOR THE PES DICT OBJECTS CONTAINING INFO FOR THE REACTIONS ON PES
def build_pes_idx_dct(pes_dct):
    """ build a dct relating index to formulat
    """
    idx_dct = {}
    form_dct = {}
    for pes_idx, formula in enumerate(pes_dct):
        idx_dct[pes_idx+1] = formula
        form_dct[formula] = pes_idx+1

    return idx_dct, form_dct


def reduce_pes_dct_to_user_inp(pes_dct, pesnums):
    """ get a pes dictionary containing only the PESs the user is running
    """
    run_pes_dct = {}
    for pes_idx, formula in enumerate(pes_dct):
        if pes_idx+1 in pesnums:
            run_pes_dct[formula] = pes_dct[formula]
    return run_pes_dct


def print_pes_channels(pes_dct):
    """ Print the PES
    """

    print('\n  Sorted Mechanism read from file:')
    for pes_idx, formula in enumerate(pes_dct):
        print('! PES:', pes_idx+1, formula)
        pes_rxn_name_lst = pes_dct[formula]['rxn_name_lst']
        pes_rct_names_lst = pes_dct[formula]['rct_names_lst']
        pes_prd_names_lst = pes_dct[formula]['prd_names_lst']
        for chn_idx, _ in enumerate(pes_rxn_name_lst):
            # print('      Channel {}: {} = {}'.format(
            #     chn_idx+1,
            #     ' + '.join(pes_rct_names_lst[chn_idx]),
            #     ' + '.join(pes_prd_names_lst[chn_idx])))
            print('  {} = {}   1.0 0.0 0.0'.format(
                ' + '.join(pes_rct_names_lst[chn_idx]),
                ' + '.join(pes_prd_names_lst[chn_idx])))

    # for (formula, pes_idx, sub_pes_idx), rxn_lst in pes_dct.items():

    #     # Print PES form and SUB PES Channels
    #     print('\nPES {}: {}, SUB PES {}'.format(
    #         pes_idx, formula, sub_pes_idx))
    #     for rxn in rxn_lst:
    #         print('  Channel {}: {} = {}'.format(
    #             rxn['chn_idx'],
    #             '+'.join(rxn['reacs']),
    #             '+'.join(rxn['prods'])))


def pes_dct_w_rxn_lsts(pes_dct, idx_dct, form_dct,
                       conn_chnls_dct, run_obj_dct):
    """ Form a new PES dictionary with the rxn_lst formatted to work
        with the drivers currently
    """
    run_pes_dct = {}
    for formula in pes_dct:

        # Set correct pes index based on the formula
        pes_idx = form_dct[formula]

        # Build the names list
        pes_rct_names_lst = pes_dct[formula]['rct_names_lst']
        pes_prd_names_lst = pes_dct[formula]['prd_names_lst']
        pes_rxn_name_lst = pes_dct[formula]['rxn_name_lst']

        # Get a list of the idxs corresponding to which channels to run
        run_chnls = []
        print('run dct', run_obj_dct)
        for pes_chn_pair in run_obj_dct:
            pes_num, chn_num = pes_chn_pair
            if idx_dct[pes_num] == formula:
                run_chnls.append(chn_num)

        # Select names from the names list corresponding to chnls to run
        # conn_chnls_dct[formula] = {sub_pes_idx: [channel_idxs]}
        for sub_pes_idx, sub_chnl_idxs in conn_chnls_dct[formula].items():
            rct_names_lst = []
            prd_names_lst = []
            rxn_name_lst = []
            rxn_model_lst = []
            rxn_chn_idxs = []
            for chn_idx in run_chnls:
                if chn_idx-1 in sub_chnl_idxs:
                    rct_names_lst.append(pes_rct_names_lst[chn_idx-1])
                    prd_names_lst.append(pes_prd_names_lst[chn_idx-1])
                    rxn_name_lst.append(pes_rxn_name_lst[chn_idx-1])
                    rxn_model_lst.append(run_obj_dct[(pes_idx, chn_idx)])
                    print('chn_idx', chn_idx)
                    rxn_chn_idxs.append(chn_idx)

            # Form reaction list (is empty if no chnls requested on sub pes)
            rxn_lst = format_run_rxn_lst(
                rct_names_lst, prd_names_lst, rxn_model_lst, rxn_chn_idxs)

            # Add the rxn lst to the pes dictionary if there is anythin
            if rxn_lst:
                run_pes_dct[(formula, pes_idx, sub_pes_idx+1)] = rxn_lst

    return run_pes_dct


def format_run_rxn_lst(rct_names_lst, prd_names_lst,
                       rxn_model_lst, rxn_chn_idxs):
    """ Get the lst of reactions to be run
    """

    # Get a list of all the species in the pes
    spc_queue = []
    for idx, _ in enumerate(rct_names_lst):
        rxn_spc = list(rct_names_lst[idx])
        rxn_spc.extend(list(prd_names_lst[idx]))
        for spc in rxn_spc:
            if spc not in spc_queue:
                spc_queue.append(spc)

    # Now loop over all the reactions to build rxn_lst
    run_lst = []
    for idx, _ in enumerate(rct_names_lst):
        spc_queue = []
        rxn_spc = list(rct_names_lst[idx])
        rxn_spc.extend(list(prd_names_lst[idx]))
        for spc in rxn_spc:
            if spc not in spc_queue:
                spc_queue.append(spc)
        run_lst.append(
            {'species': spc_queue,
             'reacs': list(rct_names_lst[idx]),
             'prods': list(prd_names_lst[idx]),
             'model': rxn_model_lst[idx],
             'chn_idx': rxn_chn_idxs[idx],
             'dummy': []}
        )

    return run_lst
