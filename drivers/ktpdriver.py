""" driver for rate constant evaluations
"""

import os
from routines.pf.rates.fit import fit_rates
from routines.pf.rates import rates as messrates
from lib.runner import rates as raterunner
from lib import printmsg


def run(pes_formula,
        spc_dct,
        thy_dct,
        rxn_lst,
        model_dct,
        run_inp_dct,
        write_messrate=True,
        run_messrate=True,
        run_fits=True):
    """ main driver for generation of full set of rate constants on a single PES
    """
    # Print the header message for the driver
    printmsg.program_header('ktp')

    # Pull stuff from dcts for now
    run_prefix = run_inp_dct['run_prefix']
    save_prefix = run_inp_dct['save_prefix']
    err_thresh = 15.
    fit_method = 'arrhenius'
    troe_param_fit_lst = ['ts1', 'ts2', 'ts3', 'alpha']

    # Build the MESS label idx dictionary for the PES
    idx_dct = messrates.make_pes_idx_dct(rxn_lst, spc_dct)

    # Set path where MESS files will be written and read
    mess_path = raterunner.get_mess_path(run_prefix, pes_formula)

    # Try and read the MESS file from the filesystem first
    mess_inp_str, dat_str_lst = raterunner.read_mess_file(mess_path)
    if mess_inp_str:
        print('mess_inp_str')
        print(mess_inp_str)

    # Write the MESS file
    if write_messrate and not mess_inp_str:
        print('Gathering elstruct data to write MESS input at {}'.format(mess_path))
        # Setting some arbitrary model to fix ene trans param passing
        # Need a better system for specifying this
        test_model = rxn_lst[0]['model']  # model for ts_0
        temps = model_dct[test_model]['options']['temps']
        pressures = model_dct[test_model]['options']['pressures']
        etrans = model_dct[test_model]['etransfer']
        pst_params = model_dct[test_model]['options']['pst_params']
        multi_info = model_dct[test_model]['options']['multi_info']
        assess_pdep = model_dct[test_model]['options']['assess_pdep']
        ene_coeff = model_dct[test_model]['options']['ene_coeff']

        # Add the irc idxs for reading for temporary things
        for spc in spc_dct:
            if 'ts_' in spc:
                spc_dct[spc]['irc_idxs'] = [
                    -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]

        print('Starting mess file preparation.')
        # Write the strings for the MESS input file
        header_str, energy_trans_str = messrates.rate_headers(
            spc_dct, rxn_lst, temps, pressures, **etrans)

        # Write the MESS strings for all the PES channels
        well_str, bim_str, ts_str, dat_lst = messrates.write_channel_mess_strs(
            spc_dct, rxn_lst, pes_formula,
            multi_info, pst_params,
            save_prefix, idx_dct,
            model_dct, thy_dct)

        # Combine strings together
        mess_inp_str = '\n'.join(
            [header_str, energy_trans_str, well_str, bim_str, ts_str])
        print('mess str')
        print(mess_inp_str)

        # Build the filesystem
        if not os.path.exists(os.path.join(run_prefix, 'MESSRATE')):
            os.mkdir(os.path.join(run_prefix, 'MESSRATE'))
        if not os.path.exists(mess_path):
            os.mkdir(mess_path)

        # Write the MESS file into the filesystem
        raterunner.write_mess_file(mess_inp_str, dat_str_lst, mess_path)

    # Run mess to produce rate output
    if run_messrate:
        raterunner.run_rates(mess_path)

    # Fit rate output to modified Arrhenius forms, print in ChemKin format
    if run_fits:
        fit_rates(spc_dct, pes_formula, idx_dct,
                  model_dct, thy_dct, ene_coeff,
                  mess_path, assess_pdep, err_thresh,
                  fit_method, troe_param_fit_lst)
