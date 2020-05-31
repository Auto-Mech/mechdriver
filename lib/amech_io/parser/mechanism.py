"""
Read the mechanism file
"""

import automol
import chemkin_io
from lib.amech_io.parser import ptt


MECH_INP = 'inp/mechanism.dat'


def build_spc_queue(rxn_lst):
    """ Build spc queue from the reaction lst for the drivers
        :return spc_queue: all the species and corresponding models in rxn
        :rtype: list[(species, model),...]
    """

    if 'all' in rxn_lst:
        # First check if rxn_lst is a bunch of species
        spc_queue = rxn_lst['all']['species']
    else:
        # Build the list from expanding the reacs and prods
        spc_queue = []
        for rxn in rxn_lst:
            model = rxn['model']
            spc_queue.extend(((reac, model) for reac in rxn['reacs']))
            spc_queue.extend(((prod, model) for prod in rxn['prods']))

    return spc_queue


def parse_mechanism_file(job_path, mech_type, spc_dct, run_obj_dct,
                         sort_rxns=False):
    """ Get the reactions and species from the mechanism input
    """
    # parsing moved to the input parsing module I am writing
    mech_str = ptt.read_inp_str(job_path, MECH_INP)
    if mech_type.lower() == 'chemkin':
        formulas, rct_names, prd_names, rxn_names = _parse_chemkin(
            mech_str, spc_dct, sort_rxns)
    # elif mech_type.lower() == 'json':
    #     spc_dct, rct_names, prd_names, rxn_name, form_strs = _parse_json(
    #         mech_path, mech_file, check_stereo, rad_rad_sort)
    else:
        raise NotImplementedError

    # Build the total PES dct
    pes_dct = build_pes_dct(formulas, rct_names, prd_names, rxn_names)

    # Build an index dct relating idx to formula
    idx_dct, form_dct = build_pes_idx_dct(pes_dct)

    # Reduce the PES dct to only what the user requests
    if run_obj_dct:
        pesnums = [idx_pair[0] for idx_pair in run_obj_dct]
    else:
        pesnums = [idx_pair[0] for idx_pair in run_obj_dct]
    reduced_pes_dct = reduce_pes_dct_to_user_inp(pes_dct, pesnums)

    # Get a dct for all of the connected channels with the PESs to run
    conn_chnls_dct = determine_connected_pes_channels(reduced_pes_dct)

    # Form the pes dct that has info formatted to run
    # Get the models in here
    run_pes_dct = pes_dct_w_rxn_lsts(
        reduced_pes_dct, idx_dct, form_dct, conn_chnls_dct, run_obj_dct)

    # Print the channels for the whole mechanism file
    print_pes_channels(pes_dct)

    return run_pes_dct


def _parse_chemkin(mech_str, spc_dct, sort_rxns):
    """ parse a chemkin formatted mechanism file
    """

    # Read the reactions and participaring species from the mech file
    rxn_block_str = chemkin_io.parser.mechanism.reaction_block(mech_str)
    rxn_strs = chemkin_io.parser.reaction.data_strings(rxn_block_str)
    rct_names_lst = list(
        map(chemkin_io.parser.reaction.reactant_names, rxn_strs))
    prd_names_lst = list(
        map(chemkin_io.parser.reaction.product_names, rxn_strs))

    # Build the inchi dct
    ich_dct = {}
    for key in spc_dct.keys():
        if 'ts' not in key and 'global' not in key:
            ich_dct[key] = spc_dct[key]['ich']

    # Sort reactant and product name lists by formula to facilitate
    # multichannel, multiwell rate evaluations
    formula_str = ''
    rxn_name_lst = []
    formula_str_lst = []
    for rct_names, prd_names in zip(rct_names_lst, prd_names_lst):
        rxn_name = '='.join(['+'.join(rct_names), '+'.join(prd_names)])
        rxn_name_lst.append(rxn_name)
        rct_ichs = list(map(ich_dct.__getitem__, rct_names))
        formula_dct = ''
        for rct_ich in rct_ichs:
            formula_i_dct = automol.inchi.formula(rct_ich)
            formula_dct = automol.formula.join(formula_dct, formula_i_dct)
        formula_str = automol.formula.string2(formula_dct)
        formula_str_lst.append(formula_str)

    rxn_info_lst = list(
        zip(formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst))

    # Sort the reactions if desired
    if sort_rxns:
        rxn_info_lst.sort()
        formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst = zip(
            *rxn_info_lst)

    return formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst


# FUNTIONS FOR THE PES DICT OBJECTS CONTAINING INFO FOR THE REACTIONS ON PES
def build_pes_dct(formula_str_lst, rct_names_lst,
                  prd_names_lst, rxn_name_lst):
    """ Build a dictionary of the PESs
    """
    pes_dct = {}
    current_formula = ''
    for fidx, formula in enumerate(formula_str_lst):
        if current_formula == formula:
            pes_dct[formula]['rct_names_lst'].append(rct_names_lst[fidx])
            pes_dct[formula]['prd_names_lst'].append(prd_names_lst[fidx])
            pes_dct[formula]['rxn_name_lst'].append(rxn_name_lst[fidx])
        else:
            current_formula = formula
            pes_dct[formula] = {}
            pes_dct[formula]['rct_names_lst'] = [rct_names_lst[fidx]]
            pes_dct[formula]['prd_names_lst'] = [prd_names_lst[fidx]]
            pes_dct[formula]['rxn_name_lst'] = [rxn_name_lst[fidx]]

    return pes_dct


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
        print('    PES:', pes_idx+1, formula)
        pes_rxn_name_lst = pes_dct[formula]['rxn_name_lst']
        pes_rct_names_lst = pes_dct[formula]['rct_names_lst']
        pes_prd_names_lst = pes_dct[formula]['prd_names_lst']
        for chn_idx, _ in enumerate(pes_rxn_name_lst):
            print('      Channel {}: {} = {}'.format(
                chn_idx+1,
                ' + '.join(pes_rct_names_lst[chn_idx]),
                ' + '.join(pes_prd_names_lst[chn_idx])))


# FUNCTIONS FOR THE CHANNELS DICT OBJECTS
def determine_connected_pes_channels(pes_dct):
    """ Determine all the connected reaction channels for each PES
        Build a dictionary for each PES with lists of connected channels:
            dct[PES_FORMULA] = [ [SUB_PES_1], [SUB_PES_2], ... , [SUB_PES_N] ]
            where each SUB_PES = [n1, n2, ... , nN],
            where n1 to nN correspond to ixds for channels that are
            connected to each other
        For efficiency we only determine channels for PESs we wish to run.
    """
    conn_chn_dct = {}
    for _, formula in enumerate(pes_dct):
        # Set the names lists for the rxns and species needed below
        pes_rct_names_lst = pes_dct[formula]['rct_names_lst']
        pes_prd_names_lst = pes_dct[formula]['prd_names_lst']
        pes_rxn_name_lst = pes_dct[formula]['rxn_name_lst']

        # Split up channels into a connected sub-pes within a formula
        subpes_idx = 0
        conndct = {}
        connchnls = {}
        for chnl_idx, _ in enumerate(pes_rxn_name_lst):
            connected_to = []
            chnl_species = [list(pes_rct_names_lst[chnl_idx]),
                            list(pes_prd_names_lst[chnl_idx])]
            for conn_chnls_idx in conndct:
                for spc_pair in chnl_species:
                    if len(spc_pair) == 1:
                        if spc_pair in conndct[conn_chnls_idx]:
                            if conn_chnls_idx not in connected_to:
                                connected_to.append(conn_chnls_idx)
                        elif spc_pair[::-1] in conndct[conn_chnls_idx]:
                            if conn_chnls_idx not in connected_to:
                                connected_to.append(conn_chnls_idx)
            if not connected_to:
                conndct[subpes_idx] = chnl_species
                connchnls[subpes_idx] = [chnl_idx]
                subpes_idx += 1
            else:
                conndct[connected_to[0]].extend(chnl_species)
                connchnls[connected_to[0]].append(chnl_idx)
                if len(connected_to) > 1:
                    for cidx, cval in enumerate(connected_to):
                        if cidx > 0:
                            conn_specs = conndct.pop(cval, None)
                            conn_chnls = connchnls.pop(cval, None)
                            conndct[connected_to[0]].extend(conn_specs)
                            connchnls[connected_to[0]].extend(conn_chnls)
                for cidx in conndct:
                    conndct[cidx].sort()
                    conndct[cidx] = [
                        conndct[cidx][i] for i in
                        range(len(conndct[cidx])) if i == 0 or
                        conndct[cidx][i] != conndct[cidx][i-1]]

        # Add connected channels list to the dictionary
        conn_chn_dct[formula] = connchnls

    return conn_chn_dct


def pes_dct_w_rxn_lsts(pes_dct, idx_dct, form_dct, conn_chnls_dct, run_obj_dct):
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
             'chn_idx': rxn_chn_idxs[idx]})

    return run_lst
# def parse_json():
#     """ parse a json file mechanism file
#     """
#
#     with open(os.path.join(mech_path, mech_file)) as mechfile:
#         inp_data = json.amech_io.parser(
#             mechfile, object_pairs_hook=collections.OrderedDict)
#         mech_data = []
#     for reaction in inp_data:
#         if isinstance(reaction, dict):
#             mech_data = inp_data
#             break
#         else:
#             for entry in inp_data[reaction]:
#                 mech_data.append(entry)
#
#     # Convert essential pieces of json file to chemkin formatted data so
#     # (i) can easily remove species that don't really exist
#     # (ii) revise products of reactions for species that don't exist
#     # (iii) do the stereochemistry generation only one
#
#     formula_str = ''
#     formula_str_lst = []
#     rxn_name_lst = []
#     rct_names_lst = []
#     prd_names_lst = []
#     rct_smis_lst = []
#     rct_ichs_lst = []
#     rct_muls_lst = []
#     prd_smis_lst = []
#     prd_ichs_lst = []
#     prd_muls_lst = []
#     prd_names_lst = []
#     rxn_sens = []
#     rxn_unc = []
#     rxn_val = []
#     rxn_fam = []
#     unq_rxn_lst = []
#     fll_rxn_lst = []
#     # idxp = 0
#     for idx, reaction in enumerate(mech_data):
#         if 'Reactants' in reaction and 'Products' in reaction:
#             print(idx, reaction['name'])
#             if reaction['name'] in fll_rxn_lst:
#                 print('duplicate reaction found:', reaction['name'], idx)
#             else:
#                 unq_rxn_lst.append(reaction['name'])
#             fll_rxn_lst.append(reaction['name'])
#     print('reaction duplicate test:', len(unq_rxn_lst), len(fll_rxn_lst))
#
#     for _, reaction in enumerate(mech_data):
#         # set up reaction info
#         rct_smis = []
#         rct_ichs = []
#         rct_muls = []
#         rct_names = []
#         prd_smis = []
#         prd_ichs = []
#         prd_muls = []
#         prd_names = []
#         if 'Reactants' in reaction and 'Products' in reaction:
#             for rct in reaction['Reactants']:
#                 rct_names.append(rct['name'])
#                 rct_smis.append(rct['SMILES'][0])
#                 ich = rct['InChi']
#                 if check_stereo:
#                     if not automol.inchi.is_complete(ich):
#                         print('adding stereochemsiry for {}'.format(ich))
#                         ich = automol.inchi.add_stereo(rct['InChi'])[0]
#                 rct_ichs.append(ich)
#                 rct_muls.append(rct['multiplicity'])
#             rad_rad_reac = True
#             if len(rct_ichs) == 1:
#                 rad_rad_reac = False
#             else:
#                 if min(rct_muls) == 1:
#                     rad_rad_reac = False
#             for prd in reaction['Products']:
#                 prd_names.append(prd['name'])
#                 prd_smis.append(prd['SMILES'][0])
#                 ich = prd['InChi']
#                 if check_stereo:
#                     if not automol.inchi.is_complete(ich):
#                         print('adding stereochemsiry for {}'.format(ich))
#                         ich = automol.inchi.add_stereo(prd['InChi'])[0]
#                 prd_ichs.append(ich)
#                 prd_muls.append(prd['multiplicity'])
#             rad_rad_prod = True
#             if len(prd_ichs) == 1:
#                 rad_rad_prod = False
#             else:
#                 if min(prd_muls) == 1:
#                     rad_rad_prod = False
#
#             print('rad_rad_sort test:', rad_rad_sort)
#             if rad_rad_sort and not rad_rad_reac and not rad_rad_prod:
#                 continue
#             rct_smis_lst.append(rct_smis)
#             rct_ichs_lst.append(rct_ichs)
#             rct_muls_lst.append(rct_muls)
#             rct_names_lst.append(rct_names)
#             prd_smis_lst.append(prd_smis)
#             prd_ichs_lst.append(prd_ichs)
#             prd_muls_lst.append(prd_muls)
#             prd_names_lst.append(prd_names)
#         rxn_name_lst.append(reaction['name'])
#         if 'Sensitivity' in reaction:
#             rxn_sens.append(reaction['sensitivity'])
#         else:
#             rxn_sens.append('')
#         if 'Uncertainty' in reaction:
#             rxn_unc.append(reaction['uncertainty'])
#         else:
#             rxn_unc.append('')
#         if 'value' in reaction:
#             rxn_val.append(reaction['value'])
#         else:
#             rxn_val.append('')
#         if 'family' in reaction:
#             rxn_fam.append(reaction['family'])
#         else:
#             rxn_fam.append('')
#
#         formula_dct = ''
#         for rct_ich in rct_ichs:
#             formula_i_dct = automol.inchi.formula(rct_ich)
#             formula_dct = automol.formula.join(
#                 formula_dct, formula_i_dct)
#         formula_str = automol.formula.string2(formula_dct)
#         formula_str_lst.append(formula_str)
#
#     print('rct_muls_lst:', rct_muls_lst)
#     unq_ich_lst = []
#     unq_mul_lst = []
#     unq_smi_lst = []
#     unq_lab_lst = []
#     unq_lab_idx_lst = []
#     csv_str = 'name,smiles,mult'
#     csv_str += '\n'
#     spc_str = 'species'
#     spc_str += '\n'
#     for ichs, muls, smis in zip(rct_ichs_lst, rct_muls_lst, rct_smis_lst):
#         for ich, mul, smi in zip(ichs, muls, smis):
#             unique = True
#             for unq_ich, unq_mul in zip(unq_ich_lst, unq_mul_lst):
#                 if ich == unq_ich and mul == unq_mul:
#                     unique = False
#             if unique:
#                 unq_ich_lst.append(ich)
#                 unq_mul_lst.append(mul)
#                 unq_smi_lst.append(smi)
#
#                 formula_dct = automol.inchi.formula(ich)
#                 lab = automol.formula.string2(formula_dct)
#
#                 unq_lab_lst.append(lab)
#                 lab_idx = -1
#                 for lab_i in unq_lab_lst:
#                     if lab == lab_i:
#                         lab_idx += 1
#                 unq_lab_idx_lst.append(lab_idx)
#                 if lab_idx == 0:
#                     label = lab
#                 else:
#                     label = lab + '(' + str(lab_idx) + ')'
#                 csv_str += ','.join([label, smi, str(mul)])
#                 csv_str += '\n'
#                 spc_str += label
#                 spc_str += '\n'
#     for ichs, muls, smis in zip(prd_ichs_lst, prd_muls_lst, prd_smis_lst):
#         for ich, mul, smi in zip(ichs, muls, smis):
#             unique = True
#             for unq_ich, unq_mul in zip(unq_ich_lst, unq_mul_lst):
#                 if ich == unq_ich and mul == unq_mul:
#                     unique = False
#             if unique:
#                 unq_ich_lst.append(ich)
#                 unq_mul_lst.append(mul)
#                 unq_smi_lst.append(smi)
#
#                 formula_dct = automol.inchi.formula(ich)
#                 lab = automol.formula.string2(formula_dct)
#
#                 unq_lab_lst.append(lab)
#                 lab_idx = -1
#                 for lab_i in unq_lab_lst:
#                     if lab == lab_i:
#                         lab_idx += 1
#                 unq_lab_idx_lst.append(lab_idx)
#                 if lab_idx == 0:
#                     label = lab
#                 else:
#                     label = lab + '(' + str(lab_idx) + ')'
#                 csv_str += ','.join([label, smi, str(mul)])
#                 csv_str += '\n'
#                 spc_str += label
#                 spc_str += '\n'
#
#     spc_str += 'END'
#     spc_str += '\n'
#     spc_str += '\n'
#
#     sort_smiles_path = os.path.join(mech_path, 'smiles_sort.csv')
#     with open(sort_smiles_path, 'w') as sorted_csv_file:
#         sorted_csv_file.write(csv_str)
#
#     rxn_info_lst = list(zip(
#        formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst, rxn_sens,
#         rxn_unc, rxn_val, rxn_fam, rct_smis_lst, rct_ichs_lst, rct_muls_lst,
#         prd_smis_lst, prd_ichs_lst, prd_muls_lst))
#     rxn_info_lst = sorted(rxn_info_lst, key=lambda x: (x[0]))
#     old_formula = rxn_info_lst[0][0]
#     sens = rxn_info_lst[0][4]
#     ordered_formula = []
#     ordered_sens = []
#     for entry in rxn_info_lst:
#         formula = entry[0]
#         if formula == old_formula:
#             sens = max(sens, entry[4])
#         else:
#             ordered_sens.append(sens)
#             ordered_formula.append(old_formula)
#             sens = entry[4]
#             old_formula = formula
#     ordered_sens.append(sens)
#     ordered_formula.append(old_formula)
#     sens_dct = {}
#     for i, sens in enumerate(ordered_sens):
#         sens_dct[ordered_formula[i]] = sens
#     rxn_info_lst = sorted(
#         rxn_info_lst, key=lambda x: (sens_dct[x[0]], x[4]), reverse=True)
#
#     formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst, rxn_sens,
#     rxn_unc, rxn_val, rxn_fam, rct_smis_lst, rct_ichs_lst, rct_muls_lst,
#     prd_smis_lst, prd_ichs_lst, prd_muls_lst = zip(*rxn_info_lst)
#
#     rxn_namep_lst = []
#     rxn_namep_str = 'REACTIONS   KCAL/MOLE   MOLES'
#     rxn_namep_str += '\n'
#     for i, _ in enumerate(rxn_name_lst):
#         rxn_namep = []
#         rct_labs = []
#         rct_dat = zip(rct_smis_lst[i], rct_ichs_lst[i], rct_muls_lst[i])
#         spc_dat = zip(unq_ich_lst, unq_mul_lst, unq_lab_lst, unq_lab_idx_lst)
#         for _, rct_ich, rct_mul in rct_dat:
#             for ich, mul, lab, lab_idx in spc_dat:
#                 if rct_ich == ich and rct_mul == mul:
#                     if lab_idx == 0:
#                         rct_lab = lab
#                     else:
#                         rct_lab = lab + '(' + str(lab_idx) + ')'
#                     break
#             rct_labs.append(rct_lab)
#         rct_label = '+'.join(rct_labs)
#         prd_labs = []
#         prd_dat = zip(prd_smis_lst[i], prd_ichs_lst[i], prd_muls_lst[i])
#         spc_dat = zip(unq_ich_lst, unq_mul_lst, unq_lab_lst, unq_lab_idx_lst)
#         for _, prd_ich, prd_mul in prd_dat:
#             for ich, mul, lab, lab_idx in spc_dat:
#                 if prd_ich == ich and prd_mul == mul:
#                     if lab_idx == 0:
#                         prd_lab = lab
#                     else:
#                         prd_lab = lab + '(' + str(lab_idx) + ')'
#                     break
#             prd_labs.append(prd_lab)
#         prd_label = '+'.join(prd_labs)
#         rate_str = str('  1.e10   1.0   10000.  ! Sens = ')
#         rxn_namep = (
#             rct_label + ' <=> ' + prd_label + rate_str + str(rxn_sens[i]))
#         rxn_namep_str += rxn_namep
#         rxn_namep_str += '\n'
#         rxn_namep_lst.append(rxn_namep)
#
#     mech_str = spc_str + rxn_namep_str
#     mech_str += 'END'
#     mech_str += '\n'
#
#   with open(os.path.join(mech_path, 'mech_sort.txt'), 'w') as sort_mech_file:
#         sort_mech_file.write(mech_str)
#
#     # set up species info
#     spc_names = []
#     chg_dct = {}
#     mul_dct = {}
#     spc_dct = {}
#     tot_lst = rct_names_lst + prd_names_lst
#     for i, spc_names_lst in enumerate(tot_lst):
#         for j, spc_name in enumerate(spc_names_lst):
#             chg = 0
#             if spc_name not in spc_names:
#                 spc_names.append(spc_name)
#                 chg_dct[spc_name] = chg
#                 print('rct_muls test1:', spc_name, i, j)
#                 print('rct_muls test2:', rct_muls_lst[i])
#                 mul_dct[spc_name] = rct_muls_lst[i][j]
#                 spc_dct[spc_name] = {}
#                 spc_dct[spc_name]['chg'] = chg
#                 spc_dct[spc_name]['ich'] = rct_ichs_lst[i][j]
#                 spc_dct[spc_name]['mul'] = rct_muls_lst[i][j]
#     for i, spc_names_lst in enumerate(prd_names_lst):
#         for j, spc_name in enumerate(spc_names_lst):
#             chg = 0
#             if spc_name not in spc_names:
#                 spc_names.append(spc_name)
#                 chg_dct[spc_name] = chg
#                 mul_dct[spc_name] = prd_muls_lst[i][j]
#                 spc_dct[spc_name] = {}
#                 spc_dct[spc_name]['chg'] = chg
#                 spc_dct[spc_name]['ich'] = prd_ichs_lst[i][j]
#                 spc_dct[spc_name]['mul'] = prd_muls_lst[i][j]
#     rxn_info_lst = list(
#         zip(formula_str_lst, rct_names_lst, prd_names_lst, rxn_name_lst))
#
#   return spc_dct, rct_names_lst, prd_names_lst, rxn_name_lst, formula_str_lst
