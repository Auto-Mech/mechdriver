""" driver for rate constant evaluations
"""

import os
from mess_io.writer import rxnchan_header_str
from routines.pf import ktp as ktp_routines
from routines.pf.runner import ktp as ktp_runner
from lib import filesys
from lib.amech_io import parser
from lib.amech_io import cleaner


def run(pes_formula, pes_idx,
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
    # Fix this to read ene model
    print('\nIdentifying reaction classes for transition states...')
    ts_dct = {}
    for rxn in rxn_lst:
        tsname = 'ts_{:g}_{:g}'.format(pes_idx, rxn['chn_idx'])
        spc_model = rxn['model'][1]
        ene_model = spc_model_dct[spc_model]['es']['ene']
        geo_model = spc_model_dct[spc_model]['es']['geo']
        es_info = parser.model.set_es_model_info(
            spc_model_dct[spc_model]['es'], thy_dct)
        if not isinstance(ene_model, str):
            ene_method = ene_model[1][1]
        else:
            ene_method = ene_model
        thy_info = filesys.inf.get_es_info(ene_method, thy_dct)
        ini_thy_info = filesys.inf.get_es_info(geo_model, thy_dct)
        pf_model = parser.model.set_pf_model_info(
            spc_model_dct[spc_model]['pf'])
        ts_dct[tsname] = parser.species.build_sing_chn_sadpt_dct(
            tsname, rxn, thy_info, ini_thy_info,
            run_inp_dct, spc_dct, cla_dct)
    spc_dct = parser.species.combine_sadpt_spc_dcts(
        ts_dct, spc_dct)

    # Build the MESS label idx dictionary for the PES
    label_dct = ktp_routines.rates.make_pes_label_dct(
        rxn_lst, pes_idx, spc_dct)

    # Set path where MESS files will be written and read
    mess_path = ktp_runner.get_mess_path(run_prefix, pes_formula)

    # Try and read the MESS file from the filesystem first
    mess_inp_str, dat_lst = ktp_runner.read_mess_file(mess_path)
    # if mess_inp_str:
    #     print('mess_inp_str')
    #     print(mess_inp_str)

    # Write the MESS file
    if write_messrate:  # and not mess_inp_str:
        print(('\n\n------------------------------------------------' +
               '--------------------------------------'))
        print('\nGathering ElStruct data to write MESS',
              'input at {}'.format(mess_path))

        # Write the strings for the MESS input file
        globkey_str, energy_trans_str = ktp_routines.rates.make_header_str(
            spc_dct, rxn_lst, temps, pressures, **etransfer)

        # Write the MESS strings for all the PES channels
        well_str, bim_str, ts_str, dat_lst = ktp_routines.rates.make_pes_mess_str(
            spc_dct, rxn_lst, pes_formula, pes_idx,
            save_prefix, label_dct,
            spc_model_dct, thy_dct)

        # Combine strings together
        rchan_header_str = rxnchan_header_str()
        mess_inp_str = '\n'.join(
            [globkey_str,
             energy_trans_str,
             rchan_header_str,
             well_str,
             bim_str,
             ts_str]
        )
        mess_inp_str = cleaner.remove_trail_whitespace(mess_inp_str)

        # Build the filesystem
        if not os.path.exists(os.path.join(run_prefix, 'MESSRATE')):
            os.mkdir(os.path.join(run_prefix, 'MESSRATE'))
        if not os.path.exists(mess_path):
            os.mkdir(mess_path)

        # Write the MESS file into the filesystem
        ktp_runner.write_mess_file(mess_inp_str, dat_lst, mess_path)

    # Run mess to produce rate output
    if run_messrate:
        print(('\n\n------------------------------------------------' +
               '--------------------------------------'))
        print('\nRunning MESS for the input file at {}'.format(mess_path))
        ktp_runner.run_rates(mess_path)
        # Include the PES number

    # Fit rate output to modified Arrhenius forms, print in ChemKin format
    if run_fits:
        print(('\n\n------------------------------------------------' +
               '--------------------------------------'))
        print('\nFitting Rate Constants for PES to Functional Forms')
        # Include the PES number
        ktp_routines.fit.fit_rates(temps, pressures, tunit, punit,
                                   pes_formula, label_dct,
                                   es_info, pf_model,
                                   mess_path, fit_method, pdep_fit,
                                   arrfit_thresh)
