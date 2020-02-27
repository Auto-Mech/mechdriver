""" driver for rate constant evaluations
"""

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
        run_rates=True,
        run_fits=True):
    """ main driver for generation of full set of rate constants on a single PES
    """
    # Print the header message for the driver
    printmsg.program_header('ktp')

    # Pull stuff from dcts for now
    save_prefix = run_inp_dct['save_prefix']

    # Run the rates
    if run_rates:

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
                    -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]

        # Getting some other info to pass
        # rct_ichs = spc_dct['ts_0']['rxn_ichs'][0]
        print('Starting mess file preparation.')

        # Write the strings for the MESS input file
        header_str, energy_trans_str = messrates.rate_headers(
            spc_dct, rxn_lst, temps, pressures, **etrans)

        # Build the MESS label idx dictionary for the PES
        idx_dct = messrates.make_pes_idx_dct(rxn_lst, spc_dct)
        print('idx dct test')
        for key, val in idx_dct.items():
            print(key, val)

        # Write the MESS strings for all the PES channels
        # idx_dct = {}
        well_str, bim_str, ts_str, dat_lst = messrates.write_channel_mess_strs(
            spc_dct, rxn_lst, pes_formula,
            multi_info, pst_params,
            save_prefix, idx_dct,
            model_dct, thy_dct)

        # Run mess to produce rate output
        mess_path = raterunner.run_rates(
            header_str, energy_trans_str,
            well_str, bim_str, ts_str, dat_lst,
            spc_dct['ts_0'], spc_dct['ts_0']['rxn_fs'][3],
            model_dct, thy_dct, test_model)

    # Fit rate output to modified Arrhenius forms, print in ChemKin format
    if run_fits:
        fit_rates(spc_dct, pes_formula, idx_dct,
                  model_dct, thy_dct, ene_coeff,
                  mess_path, assess_pdep)
