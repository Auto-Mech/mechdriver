""" driver for rate constant evaluations
"""

import os
from routines.pf.ktp.fit import fit_rates
from routines.pf.ktp import rates as messrates
from lib.runner import rates as raterunner
from lib.filesystem import inf as finf
from lib.load import species as loadspc
from lib import printmsg


def run(pes_formula, pes_idx, sub_pes_idx,
        spc_dct,
        cla_dct,
        thy_dct,
        rxn_lst,
        pes_model_dct, spc_model_dct,
        run_inp_dct,
        write_messrate=True,
        run_messrate=True,
        run_fits=True):
    """ main driver for generation of full set of rate constants on a single PES
    """
    print('rxn lst', rxn_lst)

    # Print the header message for the driver
    printmsg.program_header('ktp')

    # Pull stuff from dcts for now
    run_prefix = run_inp_dct['run_prefix']
    save_prefix = run_inp_dct['save_prefix']

    # Pull PES model and pieces
    pes_model = rxn_lst[0]['model'][0]
    temps = pes_model_dct[pes_model]['rate_temps']
    pressures = pes_model_dct[pes_model]['pressures']
    etransfer = pes_model_dct[pes_model]['etransfer']
    pdep_fit = pes_model_dct[pes_model]['pdep_fit']
    tunit = pes_model_dct[pes_model]['tunit']
    punit = pes_model_dct[pes_model]['punit']
    fit_method = pes_model_dct[pes_model]['fit_method']
    arrfit_thresh = pes_model_dct[pes_model]['dbl_arrfit_thresh']

    # Get info for the transition states (want under write..)
    for rxn in rxn_lst:
        spc_model = rxn['model'][1]
        ene_model = spc_model_dct[spc_model]['es']['ene']
        geo_model = spc_model_dct[spc_model]['es']['geo']
        # Need to fix
        # thy_info = finf.get_es_info(ene_model, thy_dct)
        thy_info = finf.get_es_info(geo_model, thy_dct)
        ini_thy_info = finf.get_es_info(geo_model, thy_dct)
        ts_dct = loadspc.build_sadpt_dct(
            pes_idx, rxn_lst, thy_info, ini_thy_info,
            run_inp_dct, spc_dct, cla_dct)
        spc_dct = loadspc.combine_sadpt_spc_dcts(
            ts_dct, spc_dct)

    # Build the MESS label idx dictionary for the PES
    label_dct = messrates.make_pes_label_dct(rxn_lst, pes_idx, spc_dct)

    # Set path where MESS files will be written and read
    mess_path = raterunner.get_mess_path(run_prefix, pes_formula)

    # Try and read the MESS file from the filesystem first
    mess_inp_str, dat_lst = raterunner.read_mess_file(mess_path)
    # if mess_inp_str:
    #     print('mess_inp_str')
    #     print(mess_inp_str)

    # Write the MESS file
    if write_messrate:  # and not mess_inp_str:
        print('Gathering ElStruct data to write MESS',
              'input at {}'.format(mess_path))

        print('Starting mess file preparation.')
        # Write the strings for the MESS input file
        globkey_str, energy_trans_str = messrates.make_header_str(
            spc_dct, rxn_lst, temps, pressures, **etransfer)

        # Write the MESS strings for all the PES channels
        well_str, bim_str, ts_str, dat_lst = messrates.make_pes_mess_str(
            spc_dct, rxn_lst, pes_formula, pes_idx, 
            save_prefix, label_dct,
            spc_model_dct, thy_dct)

        # Combine strings together
        mess_inp_str = '\n'.join(
            [globkey_str, energy_trans_str, well_str, bim_str, ts_str])

        # Build the filesystem
        if not os.path.exists(os.path.join(run_prefix, 'MESSRATE')):
            os.mkdir(os.path.join(run_prefix, 'MESSRATE'))
        if not os.path.exists(mess_path):
            os.mkdir(mess_path)

        # Write the MESS file into the filesystem
        raterunner.write_mess_file(mess_inp_str, dat_lst, mess_path)

    # Run mess to produce rate output
    if run_messrate:
        raterunner.run_rates(mess_path)

    # Fit rate output to modified Arrhenius forms, print in ChemKin format
    if run_fits:
        fit_rates(temps, pressures, tunit, punit,
                  pes_formula, label_dct,
                  mess_path, fit_method, pdep_fit,
                  arrfit_thresh)
