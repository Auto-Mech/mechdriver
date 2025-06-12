""" Return the rates from a prompt dissociation process
"""
from operator import mod
import sys
import copy
import mess_io
from mechanalyzer import calculator
from mechanalyzer.parser._util import resort_ktp_labels
from mechanalyzer.parser._util import remove_fw_rxns
from mechanalyzer.parser._util import remove_rev_rxns

def multipes_prompt_dissociation_ktp_dct(list_strs_dct, model, bf_thresh):
    """ generation of rxn_ktp_dct for multipes with unknown character
        - identifies if ped or hot
        - extracts prompt rates for all
        - provides a final non-redundant ktp dct
        :param list_strs_dct: list of the dictionaries for each rxn
                containing the corresponding file strings
        :type list_strs_dct: [{'inp': str, 'ktp_out':str, ...}, {}, ...]
    """
    if isinstance(model, list):
        model = model[0]
    # SET UP A LIST OF PED FILES AND HOT FILES DEPENDING ON THE INPUT
    rxn_type_dct = {'ped': [], 'pedhot': [], 'hot': []}
    # base ktp dct: valid for any model. contains rxns on hot PESs (unmodified by prompt)
    # and PED rxns
    rxn_ktp_dct_full0 = {}
    for strs_dct in list_strs_dct:
        ped_spc, _ = mess_io.reader.ped.ped_names(strs_dct['inp'])
        hot_spc_en = mess_io.reader.hoten.get_hot_species(strs_dct['inp'])
        if ped_spc and hot_spc_en:
            rxn_type_dct['pedhot'].append(strs_dct)
        elif ped_spc and not hot_spc_en:
            rxn_type_dct['ped'].append(strs_dct)
        elif hot_spc_en and not ped_spc:
            rxn_type_dct['hot'].append(strs_dct)

        rxn_ktp_dct_full0.update(
            mess_io.reader.rates.get_rxn_ktp_dct(
                strs_dct['ktp_out'], filter_kts=True,
                filter_reaction_types=('fake', 'self',
                                       'loss', 'capture', 'reverse'),
                #   NB KEEP THE REVERSE HERE!! ONLY WANT 1 CHNL
                relabel_reactions=True
            ))
    # resort labels to adjust prod order for different PESs with potential same prods
    rxn_ktp_dct_full0 = resort_ktp_labels(
        rxn_ktp_dct_full0)

    rxn_ktp_dct_full = rxn_chains_calc(
        rxn_type_dct['ped'], rxn_type_dct['pedhot'], rxn_type_dct['hot'],
        model, bf_thresh)

    # Add in the prompt versions of the reactions
    # remove reactions possibly written in the fw or bw direction
    ks = rxn_ktp_dct_full.keys()
    rxn_ktp_dct_full0 = remove_rev_rxns(rxn_ktp_dct_full0, ks)
    rxn_ktp_dct_full0 = remove_fw_rxns(rxn_ktp_dct_full0, ks)
    # updated dct only with necessary rxns
    rxn_ktp_dct_full.update(
        rxn_ktp_dct_full0)
    # alternatively, you can update without popping but you need to update
    # the rxn_dct_red with the other one and then replace the latter

    return rxn_ktp_dct_full


def rxn_chains_calc(PES_0, PES_i, PES_end, model, bf_thresh):
    """ generate reaction chains of successive prompt/hot PESs
    """
    def listconnected(PES_str0, PESs_str):
        retlist = []
        [retlist.append(PES_str)
         for PES_str in PESs_str if isconnected(PES_str0, PES_str)]
        return retlist

    def isconnected(ped_spc_str, hot_en_spc_str):
        """ return true if ped species is hot species of the following PES """
        connected = []
        ped_spc_all, _ = mess_io.reader.ped.ped_names(ped_spc_str['inp'])
        hot_spc_en = mess_io.reader.hoten.get_hot_species(
            hot_en_spc_str['inp'])
        for ped_spc in ped_spc_all:
            frags = []
            [frags.extend(fragset.split('+')) for fragset in ped_spc]
            connected.extend(
                list(set(list(hot_spc_en.keys())).intersection(frags)))

        return len(connected) > 0

    rxn_ktp_dct_full = {}
    # SIMPLE PED/HOT PESs: COUPLE ANY N OF PED AND HOT
    # INTERMEDIATE PEDs: JUST ALLOW 1 REACTION CHAIN FOR NOW
    for pes0 in PES_0:
        rxn_ktp_dct_ped = {}  # ped dct

        hotend = listconnected(pes0, PES_end)
        if len(hotend) > 0:
            for hot_strs_dct in hotend:
                # call the usual function of the prompt dissociation
                rxn_ktp_dct, _, _ = calculator.nonboltz.prompt_dissociation_ktp_dct(
                    pes0['inp'], pes0['ktp_out'],
                    pes0['ped'], pes0['ke_out'],
                    hot_strs_dct['inp'], hot_strs_dct['log'],
                    model, bf_thresh)

                # update dct
                rxn_ktp_dct_ped.update(
                    rxn_ktp_dct)

        pedhot = listconnected(pes0, PES_i)
        if len(pedhot) == 1:
            # create chain
            pedhot_all = copy.deepcopy(PES_i)
            pedhot_i = pedhot[0]
            pedhot_all.remove(pedhot_i)

            # PED: PES0, HOT: PEDHOT_I - EXTRACT RATES AND NEW ENERGY DISTRIBUTION
            # call the usual function of the prompt dissociation
            rxn_ktp_dct, pednew_dct, ene_bw_dct = calculator.nonboltz.prompt_dissociation_ktp_dct(
                pes0['inp'], pes0['ktp_out'],
                pes0['ped'], pes0['ke_out'],
                pedhot_i['inp'], pedhot_i['log'],
                model, bf_thresh, hot_ped_str=pedhot_i['ped'], hot_ke_out_str=pedhot_i['ke_out'])
            # you need the reaction ktp dictionary as a new input
            # you need the starting energy distribution of the products of pedhot_i
            prompt_ktp_dct = copy.deepcopy(rxn_ktp_dct)
            # rxn_ktp_dct contains all reactions in the chain,
            # prompt_ktp_dct contains only the rxns to be "expanded"

            while len(hotend) == 0:  # until you find a termination to the chain

                pedhot_new = listconnected(pedhot_i, pedhot_all)
                if len(pedhot_new) == 1:
                    pedhot_i = copy.deepcopy(pedhot_new[0])
                    pedhot_all.remove(pedhot_i)

                    prompt_ktp_dct, pednew_dct, ene_bw_dct = calculator.nonboltz.prompt_chain_ktp_dct(
                        prompt_ktp_dct, pednew_dct,
                        pedhot_i['inp'], pedhot_i['log'],
                        model, bf_thresh, ene_bw_dct, hot_ped_str=pedhot_i['ped'],
                        hot_ke_out_str=pedhot_i['ke_out'])

                    # replace rxns every time
                    rxn_ktp_dct.update(prompt_ktp_dct)

                elif len(pedhot_new) > 1:
                    print('*Error: more than 1 hot chan starting from {}. Exiting')
                    sys.exit()

                # search for end of connection
                hotend = listconnected(pedhot_i, PES_end)
                if len(hotend) > 0:
                    hot_strs_dct = hotend[0]

                    prompt_ktp_dct, _, _ = calculator.nonboltz.prompt_chain_ktp_dct(
                        prompt_ktp_dct, pednew_dct,
                        hot_strs_dct['inp'], hot_strs_dct['log'],
                        model, bf_thresh, ene_bw_dct)

                    rxn_ktp_dct.update(prompt_ktp_dct)

            # update dct
            rxn_ktp_dct_ped.update(
                rxn_ktp_dct)

        elif len(pedhot) > 1:
            print('*Error: only chains with 1 connected pedhot species PES available')
            sys.exit()

        # merge the ped dct with the full one
        # before merging: resort labels
        rxn_ktp_dct_ped = resort_ktp_labels(
            rxn_ktp_dct_ped)
        rxn_ktp_dct_full = calculator.rates.merge_rxn_ktp_dcts(
            rxn_ktp_dct_full, rxn_ktp_dct_ped
        )

    return rxn_ktp_dct_full
