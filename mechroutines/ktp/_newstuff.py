""" functinos for new stuff like non-thermal and prompyt dissoc.
"""

# import mess_io
# import ioformat


def set_prod_density_param(rgts, pesgrp_num, pes_param_dct):
    """ Figure out if densities should be calculated
    """
    if pes_param_dct is not None:
        print('prod density test')
        all_peds = pes_param_dct['peds']
        pes_peds = all_peds[pesgrp_num]

        calc_dens = (False, False)
        for ped in pes_peds:
            ped_spc = ped.split('_')[1]  # Get a string with prds
            print('ped_spc', ped_spc)
            print('rgts', rgts)
            if all(rgt in ped_spc for rgt in rgts):
                calc_dens = (True, True)
                break
    else:
        calc_dens = tuple(False for _ in rgts)

    return calc_dens


def energy_dist_params(pesgrp_num, pes_param_dct, label_dct,
                       enes=(10.0, 20.0, 30.0)):
    """ set values to determine input parameters for handling
        energy distributions in MESS calculations

        maybe just call this before the writer and pass to make_pes_str
    """

    print('label dct', label_dct)

    if pes_param_dct is not None:

        # Grab the desired PED and hot enes for the PES in the group
        all_peds = pes_param_dct['peds']
        all_hot_spc = pes_param_dct['hot']
        pes_peds = all_peds[pesgrp_num]
        pes_hot_spc = all_hot_spc[pesgrp_num]

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
        if any(pes_hot_spc):
            hot_enes_dct = {label_dct[spc]: enes for spc in pes_hot_spc}
            hot_str = ' '.join(hot_enes_dct.keys())
            print(f'Species for Hot: {hot_str}')
        else:
            hot_enes_dct = None

        # Set the micro params for writing k(E)s
        # When to set this
        micro_out_params = (0.1, 320.0, 0.1)
        print(f'Ranges for k(E) calculations: {micro_out_params}')
    else:
        ped_spc_lst = None
        hot_enes_dct = None
        micro_out_params = None

    return ped_spc_lst, hot_enes_dct, micro_out_params


def set_hot_energies(pesgrp_num, reacs, prods,
                     chnl_enes, pes_param_dct,
                     ene_range=(10.,)):
    """ Determine what hot energies should be for the requested
        species.
    """

    all_hot_spc = pes_param_dct['hot']
    pes_hot_spc = all_hot_spc[pesgrp_num]

    for spc in pes_hot_spc:
        if spc in reacs:
            ene = chnl_enes['reacs']
        if spc in prods:
            ene = chnl_enes['prods']

    enes = (ene,) + ene_range

    return enes


# def ped_therm_range():
#     """ Use the PED to determine the energy range for the hot energy
#         values
#
#         I think this will work since we loop over PESs of the PES
#         group to do write->run->read
#
#         However, we only really want to read the final prompt rates.
#         Still need to find a way to rewrite the reactions...
#     """
#     return NotImplementedError
# def prompt_process(mess_path, pes_param_dct):
#     """ Read the MESS strings and call the prompt dissoc rate code
#     """
#
#     # Get prompt dissociation parameters
#     rad_name = pes_param_dct['rad_name']
#     modeltype = pes_param_dct['modeltype']
#     bf_thresh = pes_param_dct['bf_threshold']
#
#     # Read the MESS input and output strings
#     inp_str = ioformat.read_file(mess_path, 'mess.inp')
#     ktp_out_str = ioformat.read_file(mess_path, 'mess.out')
#     ke_out_str = ioformat.read_file(mess_path, 'ke.out')
#     log_str = ioformat.read_file(mess_path, 'mess.log')
#
#     ktp_dct = prompt_dissoctiation_ktp_dct(
#         rad_name, modeltype,
#         inp_str, ktp_out_str, ke_out_str, log_str,
#         bf_threshold=bf_thresh)
#
#     return ktp_dct
# def prompt_dissoctiation_ktp_dct(rad_name, modeltype,
#                                  inp_str, ktp_out_str, ke_out_str, log_str,
#                                  bf_threshold=0.1):
#     """ Evaluate the prompt dissociation k(T,P) values.
#
#         Reads the k(T,P) values from the MESSRATE output, as well as the
#         PEDS and hot energies to modify the values accounting for prompt
#         dissociation.
#
#         :param temps: temperatures to evaluate prompt dissociation rates
#         :type temps: numpy.array
#         :param rad_name: mechanism name of radical that is dissociated
#         :type rad_name: str
#         :param modeltype: statistical energy
#             distribution model for hot radical
#         :type modeltype: str
#         :param inp_str: MESSRATE input string
#         :type inp_str: str
#         :param ktp_out_str: MESSRATE output string with k(T,P) values
#         :type ktp_out_str: str
#         :param ke_out_str: MESSRATE output string with k(E) values
#         :type ke_out_str: str
#         :param log_str: MESSRATE output .log string
#         :type log_str: str
#         :param bf_threshold: % branching fraction to include species in prods
#         :type bf_threshold: float
#         :rtype: dict[]
#     """
#     # Get the PED info
#     pedspecies, pedoutput = mess_io.reader.ped.ped_names(inp_str)
#     reacs, prods = pedspecies
#
#     # Read the rates and energies
#     ktp_dct = mess_io.reader.rates.ktp_dct(ktp_out_str, reacs, prods)
#     enes_sp, enes_ts = mess_io.reader.rates.energies(ktp_out_str)
#     _, ene_bw = mess_io.reader.rates.barriers(enes_ts, enes_sp, reacs, prods)
#
#     # 2. READ THE HOTENERGIES OUTPUT
#     spc_blocks = mess_io.reader.get_species(inp_str)
#     hot_temps,_= mess_io.reader.rates.temperatures(inp_str, mess_file='inp')
#     hot_pressures,_= mess_io.reader.rates.pressures(inp_str, mess_file='inp')
#     # drop the last element in the pressure list ('high')
#     hot_pressures = hot_pressures[:-1]
#     hotspecies = mess_io.reader.hotenergies.get_hot_names(inp_str)
#     hoten_dct = mess_io.reader.hotenergies.extract_hot_branching(
#         log_str, hotspecies,
#         list(spc_blocks.keys()), hot_temps, hot_pressures)
#
#     # 1. READ PEDOUTPUT file and reconstruct the energy distribution
#     pedspecies, pedoutput = mess_io.reader.ped.ped_names(inp_str)
#     reacs, prods = pedspecies
#     ped_df = mess_io.reader.ped.get_ped(ktp_out_str, pedspecies, enes_sp)
#     dof_info = mess_io.reader.ped.ped_dof_MW(spc_blocks[prods])
#
#     # 1b. READ ke_ped.out file and extract energy density of each fragment
#     dos_df = mess_io.reader.rates.dos_rovib(ke_out_str)
#
#     ped_df_prod1 = mess_io.reader.ped.ped_prod1(
#         ped_df, rad_name, modeltype,
#         dos_df=dos_df, dof_info=dof_info, ene_bw=ene_bw)
#
#     # 3. DERIVE T,P-PROD BRANCHING FRACTIONS; decide which species to keep
#     bf_tp_df = mess_io.reader.bf.bf_tp_df_full(
#         rad_name, ped_df_prod1, hoten_dct)
#     bf_tp_dct = mess_io.reader.bf.bf_tp_dct_filter(
#         bf_tp_df, bf_threshold, modeltype, T_all=temps)
#
#     # 4. DO ARRHENIUS FITS FOR THE SELECTED BFs
#     rxn_ktp_dct = mess_io.reader.bf.merge_bf_rates(bf_tp_dct, ktp_dct)
#
#     # rename the ktp dictionary with appropriate reaction names
#     rxn_ktp_dct = rename_ktp_dct(rxn_ktp_dct, pedspecies, label_dct)
#     # print(rxn_ktp_dct)
#
#     fitted_dct = fit_ktp_dct(rxn_ktp_dct, CWD)
#
#     return fitted_ktp_dct
# def _hot_enes_(hot_spc, all_chn_enes):
#     """ Determine the hot energies using the relative energies of PES
#     """
#     return None
# def _map_hot_names(hot_lst, label_dct):
#     """ Change the names of the hot species to the mess labels
#     """
#     return None
