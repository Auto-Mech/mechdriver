""" functinos for new stuff like non-thermal and prompyt dissoc.
"""

# import mess_io
# import ioformat


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
                 ene_range=(10.,)):
    """ Determine what hot energies should be for the requested
        species.

        Returns a dictionary where keys are the the mechanism names
        for the side of the reaction the hot spc appears and the values
        are the energies to set in the mechanism file. {side: ene_lst}
    """

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
                hot_enes_dct[side] = (ene,) + ene_range

        if not hot_enes_dct:
            hot_enes_dct = None
    else:
        hot_enes_dct = None

    return hot_enes_dct


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
 def prompt_process(mess_path, pes_param_dct):
     """ Read the MESS strings and call the prompt dissoc rate code
     """

     # Get prompt dissociation parameters
     rad_name = pes_param_dct['rad_name']
     modeltype = pes_param_dct['modeltype']
     bf_thresh = pes_param_dct['bf_threshold']

     # Read the MESS input and output strings
     inp_str = ioformat.read_file(mess_path, 'mess.inp')
     ktp_out_str = ioformat.read_file(mess_path, 'mess.out')
     ke_out_str = ioformat.read_file(mess_path, 'ke.out')
     log_str = ioformat.read_file(mess_path, 'mess.log')

     ktp_dct = prompt_dissoctiation_ktp_dct(
         rad_name, modeltype,
         inp_str, ktp_out_str, ke_out_str, log_str,
         bf_threshold=bf_thresh)

     return ktp_dct


def prompt_dissoctiation_ktp_dct(rad_name, modeltype,
                                 inp_str, ktp_out_str, ke_out_str, log_str,
                                 bf_threshold=0.1):
    """ Evaluate the prompt dissociation k(T,P) values.

        Reads the k(T,P) values from the MESSRATE output, as well as the
        PEDS and hot energies to modify the values accounting for prompt
        dissociation.

        :param temps: temperatures to evaluate prompt dissociation rates
        :type temps: numpy.array
        :param rad_name: mechanism name of radical that is dissociated
        :type rad_name: str
        :param modeltype: statistical energy
            distribution model for hot radical
        :type modeltype: str
        :param inp_str: MESSRATE input string
        :type inp_str: str
        :param ktp_out_str: MESSRATE output string with k(T,P) values
        :type ktp_out_str: str
        :param ke_out_str: MESSRATE output string with k(E) values
        :type ke_out_str: str
        :param log_str: MESSRATE output .log string
        :type log_str: str
        :param bf_threshold: % branching fraction to include species in prods
        :type bf_threshold: float
        :rtype: dict[]
    """

    # 0. EXTRACT INPUT INFORMATION from me_ped.inp
    me_ped_inp = read_file(os.path.join(CWD, OPTS['pedinput']))
    me_ped_inp = remove_comment_lines(me_ped_inp, delim_pattern=app.escape('!'))
    me_ped_inp = remove_comment_lines(me_ped_inp, delim_pattern=app.escape('#'))
    me_ped_out = read_file(os.path.join(CWD, OPTS['pedoutput']))
    species_blocks_ped = mess_io.reader.get_species(me_ped_inp)
    T_lst, _ = mess_io.reader.rates.temperatures(me_ped_inp, mess_file='inp')
    P_lst, _ = mess_io.reader.rates.pressures(me_ped_inp, mess_file='inp')
    P_lst = P_lst[:-1]  # drop the last element in the pressure list ('high')
    pedspecies, pedoutput = mess_io.reader.ped.ped_names(me_ped_inp)
    energy_dct, _, conn_lst_dct, _ = mess_io.reader.pes(me_ped_inp)
    
    # 1. INFO FROM rate_ped.out and ke_ped.out: rate dct, energy barriers, dofs, fragment names
    ktp_dct = {}
    E_BW_dct = {}
    dof_dct = {}
    fragments_dct = {}
    # get the rates for all set of pedspecies
    for species in pedspecies:
        reacs, prods = species
        label = '->'.join(species)
        ktp_dct[label] = mess_io.reader.rates.ktp_dct(
            me_ped_out, reacs, prods)
        # find the corresponding energy barrier
        barrier_label = mess_io.reader.find_barrier(conn_lst_dct, reacs, prods)
        try:
            E_BW_dct[label] = energy_dct[barrier_label]-energy_dct[prods]
        except KeyError:
            E_BW_dct[label] = energy_dct[reacs]-energy_dct[prods]
        # derive dofs involved
        dof_info = mechanalyzer.calculator.statmodels.get_dof_info(species_blocks_ped[prods], ask_for_ts=True)
        dof_dct[label] = dof_info
        fragments_dct[label] = mess_io.reader.dct_species_fragments(species_blocks_ped)[prods]
    
    # 2. read PED
    pedoutput_str = read_file(os.path.join(CWD, pedoutput))
    ped_dct = mess_io.reader.ped.get_ped(pedoutput_str, pedspecies, energy_dct)
    
    # 3. READ THE ke_ped.out file and extract the energy density of each fragment
    ke_ped_out = read_file(os.path.join(CWD, OPTS['pedoutputmicro']))
    dos_df = mess_io.reader.rates.dos_rovib(ke_ped_out)
    
    # 4. READ THE HOTENERGIES OUTPUT
    hot_inp = read_file(os.path.join(CWD, OPTS['hotinput']))
    hot_out = read_file(os.path.join(CWD, OPTS['hotoutput']))
    species_blocks_hoten = mess_io.reader.get_species(hot_inp)
    hot_frag_dct = mess_io.reader.dct_species_fragments(species_blocks_hoten)
    T_lst_hot, _ = mess_io.reader.rates.temperatures(hot_inp, mess_file='inp')
    P_lst_hot, _ = mess_io.reader.rates.pressures(hot_inp, mess_file='inp')
    P_lst_hot = P_lst_hot[:-1] #drop last value of pressure
    hotspecies = mess_io.reader.hoten.get_hot_names(hot_inp)
    hoten_dct = mess_io.reader.hoten.extract_hot_branching(
        hot_out, hotspecies, list(species_blocks_hoten.keys()), T_lst_hot, P_lst_hot)
    
    # DERIVE BF AND RATES
    rxns = {}
    for species in pedspecies:
        label = '->'.join(species)
        ped_df = ped_dct[label]
        E_BW = E_BW_dct[label]
        # select the fregment of which you want the PED: it is the one in common with hotspecies
        fragments = fragments_dct[label]
        try:
            frag1 = list(set(hotspecies).intersection(fragments))[0]
            fragments.remove(frag1)
            frag2 = fragments[0]
        except IndexError:
            print('no superposition between PED fragments and hot fragments - exiting now \n')
            sys.exit()
        # DERIVE PED OF THE HOT FRAGMENT
        ped_df_frag1_dct = mechanalyzer.builder.ped.ped_frag1(
            ped_df, frag1, frag2, modeltype_list, dos_df=dos_df, dof_info=dof_dct[label], E_BW=E_BW)
    
        # JOIN PED AND HOTEN -> DERIVE PRODUCTS BF
        bf_tp_dct = mechanalyzer.builder.bf.bf_tp_dct(
            modeltype_list, ped_df_frag1_dct, hoten_dct[frag1], bf_threshold, savefile=True)
    
        # NEW KTP DICTIONARY
        frag_reacs = mess_io.reader.dct_species_fragments(species_blocks_ped)[species[0]]
        rxn_ktp_dct = mechanalyzer.builder.bf.merge_bf_ktp(bf_tp_dct, ktp_dct[label], frag_reacs, frag1, frag2, hot_frag_dct)
        rxns[label] = rxn_ktp_dct

     return fitted_ktp_dct
