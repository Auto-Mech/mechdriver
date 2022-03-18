""" functinos for new stuff like non-thermal and prompyt dissoc.
"""

import numpy
import mess_io
import mechanalyzer
from mechlib.amech_io import reader


def set_prod_density_param(rgts, pesgrp_num, pes_param_dct):
    """ Figure out if densities should be calculated
    """
    if pes_param_dct is not None:
        all_peds = pes_param_dct['peds']
        pes_peds = all_peds[pesgrp_num]

        calc_dens = (False, False)
        for ped in pes_peds:
            ped_spc = ped.split('_')[1]  # Get a string with prds
            if all(rgt in ped_spc for rgt in rgts):
                calc_dens = (True, True)
                break
    else:
        calc_dens = tuple(False for _ in rgts)

    return calc_dens


def energy_dist_params(pesgrp_num, pes_param_dct, hot_enes_dct):
    """ set values to determine input parameters for handling
        energy distributions in MESS calculations

        maybe just call this before the writer and pass to make_pes_str
    """

    if pes_param_dct is not None:

        # Grab the desired PED and hot enes for the PES in the group
        # Currenly just printing them, may need to move
        all_peds = pes_param_dct['peds']
        pes_peds = all_peds[pesgrp_num]

        # Set the PEDs
        if any(pes_peds):
            ped_spc_lst = pes_peds
            ped_str = ' '.join(ped_spc_lst)
            print(f'Species for PED: {ped_str}')
        else:
            ped_spc_lst = None

        # Set the Hot Energies section
        if hot_enes_dct is not None:
            hot_str = ' '.join(hot_enes_dct.keys())
            print(f'Species for Hot: {hot_str}')

        # Set the micro params for writing k(E)s (When to set this?)
        micro_out_params = (0.1, 320.0, 0.1)
        print(f'Ranges for k(E) calculations: {micro_out_params}')
    else:
        ped_spc_lst, micro_out_params = None, None

    return ped_spc_lst, micro_out_params


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
                                pes_mod_dct, pes_mod,
                                tsk_key_dct):
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
            pes_param_dct, rate_strs_dct, mess_paths_dct,
            pes_mod_dct[pes_mod]['rate_temps'],
            pes_mod_dct[pes_mod]['pressures']
        )

    return rxn_ktp_dct


def _single_pes_ktp_dct(pes_grp_rlst,
                        tsk_key_dct, rate_strs_dct, mess_paths_dct,
                        temps, pressures):
    """ Read the rates from a single PES
    """

    # Read options
    mess_version = tsk_key_dct['mess_version']
    use_well_extension = tsk_key_dct['well_extension']

    # Read the appropriate file based on desired versions
    pes_inf = tuple(pes_grp_rlst.keys())[0]

    if mess_version == 'v1' and use_well_extension:
        typ = f'wext-{mess_version}'
        mess_path = mess_paths_dct[pes_inf][typ]
        mess_str = rate_strs_dct[pes_inf][typ]['ktp_out']
    elif mess_version == 'v1' and not use_well_extension:
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
    else:
        rxn_ktp_dct = None
        print(f'No MESS output found at {mess_path}')

    return rxn_ktp_dct


def _prompt_dissociation_ktp_dct(pes_grp_rlst,
                                 pes_param_dct, rate_strs_dct, mess_paths_dct,
                                 temps, pressures):
    """ Evaluate the prompt dissociation k(T,P) values.

        Reads the k(T,P) values from the MESSRATE output, as well as the
        PEDS and hot energies to modify the values accounting for prompt
        dissociation.

        :param bf_threshold: % branching fraction to include species in prods
        :type bf_threshold: float
        :rtype: dict[]
    """

    # Get the PES info objects for the PED and Hot surface
    ped_pes_inf = tuple(pes_grp_rlst.keys())[0]
    hot_pes_inf = tuple(pes_grp_rlst.keys())[1]

    # Obtain the strings that are needed
    ped_mess_path = mess_paths_dct[ped_pes_inf]['base-v1']
    hot_mess_path = mess_paths_dct[hot_pes_inf]['base-v1']
    ped_strs_dct = rate_strs_dct[ped_pes_inf]['base-v1']
    hot_strs_dct = rate_strs_dct[hot_pes_inf]['base-v1']

    print('Fitting rates from\n'
          f'  - PED: {ped_mess_path}\n'
          f'  - HOT: {hot_mess_path}')

    return mechanalyzer.calculator.prompt_dissociation_ktp_dct(
        ped_strs_dct['inp'], ped_strs_dct['ktp_out'],
        ped_strs_dct['ped'], ped_strs_dct['ke_out'],
        hot_strs_dct['inp'], hot_strs_dct['ktp_out'], hot_strs_dct['log'],
        pes_param_dct['modeltype'], pes_param_dct['bf_threshold'])
# add to function as None options: temps, pressures)
