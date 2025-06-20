""" functions for new stuff like non-thermal and prompyt dissoc.
"""

import numpy
import mess_io
import mechanalyzer
from mechlib.amech_io import reader
from mechroutines.models.typ import is_abstraction_pes


def set_prod_density_param(rgts, pesgrp_num, pes_param_dct):
    """ Figure out if densities should be calculated
    """
    if pes_param_dct is not None:
        all_peds = pes_param_dct['peds']
        pes_peds = all_peds[pesgrp_num]

        calc_dens = (False, False)
        for ped in pes_peds:
            ped_spc = ped.split('=')[1]  # Get a string with prds
            if all(rgt in ped_spc for rgt in rgts):
                calc_dens = (True, True)
                break
    else:
        calc_dens = tuple(False for _ in rgts)

    return calc_dens

def relabel_ped_spc_lst(ped_spc_lst):
    """ check that all ped spc lists are of type 'A+B=C+D'.
        if A=C+D or A+B=C is found: rewrite as C+D=C+D, A+B=A+B
        reason: mess does not accept PEDs with unimol
    """
    ped_spc_lst_new = []
    for ped in ped_spc_lst:
        reacs, prods = ped.split('=')
        if '+' not in reacs:
            ped_spc_lst_new.append('='.join([prods, prods]))
        elif '+' not in prods:
            ped_spc_lst_new.append('='.join([reacs, reacs]))
        else:
            ped_spc_lst_new.append(ped)

    return ped_spc_lst_new

def energy_dist_params(pesgrp_num, pes_param_dct, hot_enes_dct, rxn_chan_str):
    """ set values to determine input parameters for handling
        energy distributions in MESS calculations

        maybe just call this before the writer and pass to make_pes_str
    """

    micro_limit = 320.  # default

    def get_ped_ene_info(pes_peds, rxn_chan_str):
        """ info for PED species - get approximate energy limit
        """
        spc_blocks_ped = mess_io.reader.get_species(rxn_chan_str)
        energy_dct, _, _, _ = mess_io.reader.pes(rxn_chan_str)
        max_ene = []
        max_ene_ped = []
        print(rxn_chan_str)
        print('peds', pes_peds)
        for ped in pes_peds:
            reacs, prods = ped.split('=')
            print('ene dct test')
            print(energy_dct)
            try:
                ene_bw = energy_dct[reacs] - energy_dct[prods]
            except KeyError as keyerr:
                print('species not in pes:', keyerr.args[0])
                print('probably unstable- setting default max_en as 300')
                max_ene.append(300.1)
                max_ene_ped.append(300.1)
                continue
            dof_info = mechanalyzer.calculator.spinfo_frommess.get_info(
                spc_blocks_ped[prods])
            max_ene_ped.append(mechanalyzer.calculator.ene_util.max_en_auto(
                dof_info['n_atoms']['TS'], ene_bw, ref_ene=energy_dct[prods]))
            max_ene.append(mechanalyzer.calculator.ene_util.max_en_auto(
                dof_info['n_atoms']['TS'], ene_bw))
            
        micro_limit = max(max_ene_ped)
        
        return max_ene, micro_limit

    if pes_param_dct is not None:

        # Grab the desired PED and hot enes for the PES in the group
        # Currently just printing them, may need to move
        all_peds = pes_param_dct['peds']
        pes_peds = all_peds[pesgrp_num]

        # Set the PEDs
        if any(pes_peds):
            ped_spc_lst = pes_peds
            ped_str = ' '.join(ped_spc_lst)
            print(f'Species for PED: {ped_str}')
            max_ene, micro_limit = get_ped_ene_info(pes_peds, rxn_chan_str)
            pes_param_dct['en_limit'][pesgrp_num] = max_ene
        else:
            ped_spc_lst = None

        # Set the Hot Energies section
        if hot_enes_dct is not None:
            hot_str = ' '.join(hot_enes_dct.keys())
            print(f'Species for Hot: {hot_str}')
            # set limit based on hoten
            micro_limit = max([max(hoten) for hoten in hot_enes_dct.values()])
        # Set the micro params for writing k(E)s (When to set this?)

        micro_out_params = (0.1, micro_limit, 0.1)
        print(f'Ranges for k(E) calculations: {micro_out_params}')
    else:
        ped_spc_lst, micro_out_params = None, None
    # relabel ped_spc_lst if needed
    if ped_spc_lst is not None:
        ped_spc_lst = relabel_ped_spc_lst(ped_spc_lst)
        
    return ped_spc_lst, micro_out_params, pes_param_dct

def set_hot_enes(hot_enes_dct, pesgrp_num, reacs, prods,
                 chnl_enes, pes_param_dct):
    """ Determine what hot energies should be for the requested
        species.

        Returns a dictionary where keys are the the mechanism names
        for the side of the reaction the hot spc appears and the values
        are the energies to set in the mechanism file. {side: ene_lst}
    """
    def set_ene_max(spc, pes_param_dct):
        """ Determine max en to write as hoten max """
        max_ene_lst = []
        for grp_num, ped_max in enumerate(pes_param_dct['en_limit']):
            for ped_i, en in enumerate(ped_max):
                if spc in pes_param_dct['peds'][grp_num][ped_i].split('=')[1]:
                    max_ene_lst.append(en)

        return max(max_ene_lst)
    # default ene_range

    if pes_param_dct is not None:
        all_hot_spc = pes_param_dct['hot']
        pes_hot_spc = all_hot_spc[pesgrp_num]

        for spc in pes_hot_spc:
            if spc in reacs:
                ene = chnl_enes['reacs']
                side = '+'.join(reacs)
            elif spc in prods:
                ene = chnl_enes['prods']
                side = '+'.join(prods)
            else:
                side, ene = None, None

            # update ene_max val if possible
            if side is not None:
                ene_max = set_ene_max(spc, pes_param_dct)
                ene_range = numpy.array([0.0 + ene, 1.0, ene_max + ene])
                hot_enes_dct[side] = tuple(ene_range)
                print('Setting {:.1f} as max hoten val for {} \n'.format(
                    ene_max+ene, side))

    return hot_enes_dct

def obtain_multipes_rxn_ktp_dct(pes_grp_rlst,
                                rate_paths_dct, pes_param_dct,
                                pes_mod_dct, pes_mod,
                                tsk_key_dct, spc_dct):
    """ Obtain the rate constants for all of the PESs in the group.
        Call additional
    """

    # Read the MESS input and output for all PES group members
    rate_strs_dct, mess_paths_dct = reader.mess.rate_strings(rate_paths_dct)

    # Read MESS file and get rate constants
    if len(pes_grp_rlst) == 1:
        rxn_ktp_dct = _single_pes_ktp_dct(
            pes_grp_rlst,
            tsk_key_dct, rate_strs_dct, mess_paths_dct,
            pes_mod_dct[pes_mod]['rate_temps'],
            pes_mod_dct[pes_mod]['pressures']
        )
    else:
        rxn_ktp_dct = _prompt_dissociation_ktp_dct(
            pes_grp_rlst,
            tsk_key_dct, spc_dct,
            pes_param_dct, rate_strs_dct, mess_paths_dct
        )

    return rxn_ktp_dct

def _single_pes_ktp_dct(pes_grp_rlst,
                        tsk_key_dct, rate_strs_dct, mess_paths_dct,
                        temps, pressures):
    """ Read the rates from a single PES
    """

    # Read options
    mess_version = tsk_key_dct['mess_version']
    use_well_lumping = tsk_key_dct['well_lumping']

    # Read the appropriate file based on desired versions
    pes_inf = tuple(pes_grp_rlst.keys())[0]

    if mess_version == 'v1' and use_well_lumping:
        typ = f'wext-{mess_version}'
        mess_path = mess_paths_dct[pes_inf][typ]
        mess_str = rate_strs_dct[pes_inf][typ]['ktp_out']
    elif mess_version == 'v1' and not use_well_lumping:
        typ = f'base-{mess_version}'
        mess_path = mess_paths_dct[pes_inf][typ]
        mess_str = rate_strs_dct[pes_inf][typ]['ktp_out']
    elif mess_version == 'v2':
        typ = f'base-{mess_version}'
        mess_path = mess_paths_dct[pes_inf][typ]
        mess_str = rate_strs_dct[pes_inf][typ]['ktp_out']

    # If file found, read and fit the rate constants
    if mess_str is not None:
        print('Fitting rates for single PES...')
        print(f'Fitting rates from {mess_path}')

        rxn_ktp_dct = mess_io.reader.rates.get_rxn_ktp_dct(
            mess_str,
            filter_kts=True,
            filter_reaction_types=('fake', 'self',
                                   'loss', 'capture', 'reverse'),
            tmin=min(temps),
            tmax=max(temps),
            pmin=min(pressures),
            pmax=max(pressures)
        )
        print('Reaction dict test')
        print(rxn_ktp_dct.keys())
    else:
        rxn_ktp_dct = None
        print(f'No MESS output found at {mess_path}')

    return rxn_ktp_dct

def _prompt_dissociation_ktp_dct(pes_grp_rlst,
                                 tsk_key_dct, spc_dct,
                                 pes_param_dct, rate_strs_dct, mess_paths_dct):
    """ Evaluate the prompt dissociation k(T,P) values.

        Reads the k(T,P) values from the MESSRATE output, as well as the
        PEDS and hot energies to modify the values accounting for prompt
        dissociation.

        :param bf_threshold: % branching fraction to include species in prods
        :type bf_threshold: float
        :rtype: dict[]
    """

    all_mess_paths = []
    strs_dct_lst = []
    for (pes_inf, rxn_lst) in pes_grp_rlst.items():
        _, pes_idx, _ = pes_inf
        if (
            tsk_key_dct['well_lumping'] and
            not is_abstraction_pes(spc_dct, rxn_lst, pes_idx)
        ):
            typ = 'wext-v1'
        else:
            typ = 'base-v1'
        all_mess_paths.append(mess_paths_dct[pes_inf][typ])
        strs_dct_lst.append(rate_strs_dct[pes_inf][typ])

    print('Fitting rates from')
    for path in all_mess_paths:
        print(f'{path}')

    return mechanalyzer.builder.multipes_prompt_dissociation_ktp_dct(
        strs_dct_lst,
        pes_param_dct['modeltype'], pes_param_dct['bf_threshold'])
