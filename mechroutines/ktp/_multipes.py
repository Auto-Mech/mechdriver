""" functinos for new stuff like non-thermal and prompyt dissoc.
"""

import numpy
import automol
import mess_io
import mechanalyzer
from mechlib.amech_io import reader


def set_prod_density_param(rgts, pesgrp_num, pes_param_dct):
    """ Figure out if densities should be calculated
    """
    if pes_param_dct is not None:
        # print('prod density test')
        all_peds = pes_param_dct['peds']
        pes_peds = all_peds[pesgrp_num]

        calc_dens = (False, False)
        for ped in pes_peds:
            ped_spc = ped.split('_')[1]  # Get a string with prds
            # print('ped_spc', ped_spc)
            # print('rgts', rgts)
            if all(rgt in ped_spc for rgt in rgts):
                calc_dens = (True, True)
                break
    else:
        calc_dens = tuple(False for _ in rgts)

    return calc_dens


def energy_dist_params(pesgrp_num, pes_param_dct, hot_enes_dct, label_dct):
    """ set values to determine input parameters for handling
        energy distributions in MESS calculations

        maybe just call this before the writer and pass to make_pes_str
    """

    if pes_param_dct is not None:

        # Grab the desired PED and hot enes for the PES in the group
        all_peds = pes_param_dct['peds']
        pes_peds = all_peds[pesgrp_num]

        # Set the PEDs
        if any(pes_peds):
            ped_spc_lst = tuple()
            for ped in pes_peds:
                _ped = ped.split('_')
                ped_spc_lst += (f'{label_dct[_ped[0]]}_{label_dct[_ped[1]]}',)
            ped_str = ' '.join(ped_spc_lst)
            print(f'Species for PED: {ped_str}')
        else:
            ped_spc_lst = None

        # Set the Hot Energies section
        if hot_enes_dct is not None:
            _hot_enes_dct = {label_dct[spc]: enes
                             for spc, enes in hot_enes_dct.items()}
            hot_str = ' '.join(_hot_enes_dct.keys())
            print(f'Species for Hot: {hot_str}')
        else:
            _hot_enes_dct = None

        # Set the micro params for writing k(E)s
        # When to set this
        micro_out_params = (0.1, 320.0, 0.1)
        print(f'Ranges for k(E) calculations: {micro_out_params}')
    else:
        ped_spc_lst = None
        _hot_enes_dct = None
        micro_out_params = None

    return ped_spc_lst, _hot_enes_dct, micro_out_params


def set_hot_enes(pesgrp_num, reacs, prods,
                 chnl_enes, pes_param_dct,
                 ene_range=None):
    """ Determine what hot energies should be for the requested
        species.

        Returns a dictionary where keys are the the mechanism names
        for the side of the reaction the hot spc appears and the values
        are the energies to set in the mechanism file. {side: ene_lst}
    """

    if ene_range is None:
        ene_range = numpy.arange(0.0, 226.0, 1.0).tolist()

    if pes_param_dct is not None:
        all_hot_spc = pes_param_dct['hot']
        pes_hot_spc = all_hot_spc[pesgrp_num]

        hot_enes_dct = {}
        for spc in pes_hot_spc:
            if spc in reacs:
                ene = chnl_enes['reacs']
                side = '+'.join(reacs)
            elif spc in prods:
                ene = chnl_enes['prods']
                side = '+'.join(prods)
            else:
                side, ene = None, None

            if side is not None:
                hot_enes_dct[side] = tuple(ene+x for x in ene_range)

        if not hot_enes_dct:
            hot_enes_dct = None
    else:
        hot_enes_dct = None

    return hot_enes_dct


def obtain_multipes_rxn_ktp_dct(pes_grp_rlst,
                                rate_paths_dct, pes_param_dct,
                                label_dct, pes_mod_dct, pes_mod):
    """ Obtain the rate constants for all of the PESs in the group.
        Call additional
    """

    # Read the MESS input and output for all PES group members
    rate_strs_dct, mess_paths_dct = reader.mess.rate_strings(rate_paths_dct)

    # Get an inverted label dct mess->mech
    # inv_label_dct = automol.util.dict_.invert(label_dct)

    # Read MESS file and get rate constants
    if len(pes_grp_rlst) == 1:

        pes_inf = tuple(pes_grp_rlst.keys())[0]

        rate_str_dct = rate_strs_dct[pes_inf]
        if rate_str_dct['wext']:
            mess_path = mess_paths_dct[pes_inf]['wext']
            mess_str = rate_str_dct['wext']['ktp_out']
        elif rate_str_dct['base']:
            mess_path = mess_paths_dct[pes_inf]['base']
            mess_str = rate_str_dct['base']['ktp_out']
        else:
            print('No Rates to fit')

        print('Fitting rates for single PES...')
        print(f'Fitting rates from {mess_path}')

        # rxn_ktp_dct = mess_io.reader.rates.get_rxn_ktp_dct(
        rxn_ktp_dct = mess_io.reader.new_rates.get_rxn_ktp_dct(
            mess_str,
            filter_kts=True,
            read_fake=False,
            read_self=False,
            read_rev=False,
            tmin=min(pes_mod_dct[pes_mod]['rate_temps']),
            tmax=max(pes_mod_dct[pes_mod]['rate_temps']),
            pmin=min(pes_mod_dct[pes_mod]['pressures']),
            pmax=max(pes_mod_dct[pes_mod]['pressures'])
        )
        # rxn_ktp_dct = mess_io.reader.rates.relabel(rxn_ktp_dct, inv_label_dct)
    else:
        rxn_ktp_dct = prompt_dissociation_ktp_dct(
            pes_grp_rlst,
            pes_param_dct, rate_strs_dct, mess_paths_dct,
            inv_label_dct,
            pes_mod_dct[pes_mod]['rate_temps'],
            pes_mod_dct[pes_mod]['pressures']
        )

    return rxn_ktp_dct


def prompt_dissociation_ktp_dct(pes_grp_rlst,
                                pes_param_dct, rate_strs_dct, mess_paths_dct,
                                inv_label_dct,
                                temps, pressures):
    """ Evaluate the prompt dissociation k(T,P) values.

        Reads the k(T,P) values from the MESSRATE output, as well as the
        PEDS and hot energies to modify the values accounting for prompt
        dissociation.

        :param bf_threshold: % branching fraction to include species in prods
        :type bf_threshold: float
        :rtype: dict[]
    """

    # Get prompt dissociation parameters
    # rad_name = pes_param_dct['rad_name']
    modeltype = pes_param_dct['modeltype']
    bf_thresh = pes_param_dct['bf_threshold']

    # Get the strings and paths
    ped_pes_inf = tuple(pes_grp_rlst.keys())[0]
    hot_pes_inf = tuple(pes_grp_rlst.keys())[0]

    # Obtain the strings that are needed
    ped_mess_path = mess_paths_dct[ped_pes_inf]['base']
    hot_mess_path = mess_paths_dct[hot_pes_inf]['base']
    print('Fitting rates from\n'
          f'  - PED: {ped_mess_path}\n'
          f'  - HOT: {hot_mess_path}')

    ped_strs_dct = rate_strs_dct[ped_pes_inf]['base']
    ped_inp_str, ped_out_str = ped_strs_dct['inp'], ped_strs_dct['ktp_out']
    ped_ped_str = ped_strs_dct['ped']
    ped_ke_out_str = ped_strs_dct['ke_out']

    hot_strs_dct = rate_strs_dct[hot_pes_inf]['base']
    hot_inp_str, hot_out_str = hot_strs_dct['inp'], hot_strs_dct['ktp_out']
    hot_log_str = hot_strs_dct['log']

    # 0. EXTRACT INPUT INFORMATION from me_ped.inp
    spc_blocks_ped = mess_io.reader.get_species(ped_inp_str)
    ped_spc, _ = mess_io.reader.ped.ped_names(ped_inp_str)  # can supply
    energy_dct, _, conn_lst_dct, _ = mess_io.reader.pes(ped_inp_str)

    # 0b. Read the rate constants from the two PESs
    rxn_ktp_dct = {}
    for mess_str in (ped_out_str, hot_out_str):
        rxn_ktp_dct.update(
            mess_io.reader.rates.get_rxn_ktp_dct(
                mess_str,
                filter_kts=True,
                tmin=min(temps),
                tmax=max(temps),
                pmin=min(pressures),
                pmax=max(pressures)
            )
        )

    # 1. INFO FROM rate_ped.out and ke_ped.out:
    #      rate dct, energy barriers, dofs, fragment names
    ene_bw_dct = {}
    dof_dct = {}
    fragments_dct = {}
    # get the rates for all set of pedspecies (using the old labels)
    for spc in ped_spc:
        reacs, prods = spc
        label = ((reacs,), (prods,), (None,))

        # Find the corresponding energy barrier
        barrier_label = mess_io.reader.find_barrier(conn_lst_dct, reacs, prods)
        try:
            ene_bw_dct[label] = energy_dct[barrier_label]-energy_dct[prods]
        except KeyError:
            ene_bw_dct[label] = energy_dct[reacs]-energy_dct[prods]

        # Derive dofs involved
        dof_info = mechanalyzer.calculator.statmodels.get_dof_info(
            spc_blocks_ped[prods], ask_for_ts=True)
        dof_dct[label] = dof_info
        fragments_dct[label] = mess_io.reader.dct_species_fragments(
            spc_blocks_ped)[prods]

    # 2. read PED
    ped_dct = mess_io.reader.ped.get_ped(ped_ped_str, ped_spc, energy_dct)

    # 3. READ ke_ped.out file and extract the energy density of each fragment
    dos_df = mess_io.reader.rates.dos_rovib(ped_ke_out_str)

    # 4. READ THE HOTENERGIES OUTPUT
    spc_blocks_hoten = mess_io.reader.get_species(hot_inp_str)
    hot_frag_dct = mess_io.reader.dct_species_fragments(spc_blocks_hoten)
    hot_spc = mess_io.reader.hoten.get_hot_names(hot_inp_str)  # can supply
    hoten_dct = mess_io.reader.hoten.extract_hot_branching(
        hot_log_str, hot_spc, list(spc_blocks_hoten.keys()),
        temps, pressures)

    # DERIVE BF AND RATES
    prompt_rxns = ()
    prompt_rxn_ktp_dct = {}
    for spc in ped_spc:
        reacs, prods = spc
        label = ((reacs,), (prods,), (None,))
        _ped_label = '+'.join(label[0]) + '->' + '+'.join(label[1])

        ped_df = ped_dct[_ped_label]
        ene_bw = ene_bw_dct[label]
        # select the fragment of which you want the PED:
        # it is the one in common with hotspecies
        fragments = fragments_dct[label]
        try:
            frag1 = list(set(hot_spc).intersection(fragments))[0]
            frag2 = list(set(fragments).difference((frag1,)))[0]
        except IndexError:
            print('no superposition between PED fragments and hot fragments '
                  '- exiting now \n')
        # DERIVE PED OF THE HOT FRAGMENT
        ped_df_frag1_dct = mechanalyzer.builder.ped.ped_frag1(
            ped_df, frag1, frag2, (modeltype,),
            dos_df=dos_df, dof_info=dof_dct[label], ene_bw=ene_bw)

        # JOIN PED AND HOTEN -> DERIVE PRODUCTS BF
        bf_tp_dct = mechanalyzer.builder.bf.bf_tp_dct(
            (modeltype,), ped_df_frag1_dct, hoten_dct[frag1], bf_thresh,
            savefile=False)

        # NEW KTP DICTIONARY
        frag_reacs_dct = mess_io.reader.dct_species_fragments(
            spc_blocks_ped)
        frag_reacs = frag_reacs_dct[spc[0]]

        prompt_rxn_ktp_dct.update(
            mechanalyzer.builder.bf.merge_bf_ktp(
                bf_tp_dct, rxn_ktp_dct[label],
                frag_reacs, frag1, frag2, hot_frag_dct)[modeltype]
        )
        prompt_rxns += (label,)

    # Remove the original reaction (currently in MESS labels
    for rxn in prompt_rxns:
        rxn_ktp_dct.pop(rxn)

    # Add in the prompt versions of the reactions
    rxn_ktp_dct.update(prompt_rxn_ktp_dct)

    # Remap the names of the remaining, non-prompt reactions ktp dct
    rxn_ktp_dct = mess_io.reader.rates.relabel(rxn_ktp_dct, inv_label_dct)

    return rxn_ktp_dct
