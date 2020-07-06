""" driver for rate constant evaluations
"""

import os
from routines.pf import ktp as ktproutines
from routines.pf import runner as pfrunner
from lib import filesys
from lib.amech_io import parser


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
    arrfit_thresh = (
        pes_model_dct[pes_model]['dbl_arrfit_thresh'],
        'max'
        # pes_model_dct[pes_model]['dbl_arrfit_val?']
    )

    # Fix this to read ene model
    print('\nIdentifying reaction classes for transition states...')
    ts_dct = {}
    for rxn in rxn_lst:
        tsname = 'ts_{:g}_{:g}'.format(pes_idx, rxn['chn_idx'])
        spc_model = rxn['model'][1]
        ene_model = spc_model_dct[spc_model]['es']['ene']
        geo_model = spc_model_dct[spc_model]['es']['geo']
        es_info = parser.model.pf_level_info(
            spc_model_dct[spc_model]['es'], thy_dct)
        if not isinstance(ene_model, str):
            ene_method = ene_model[1][1]
        else:
            ene_method = ene_model
        thy_info = filesys.inf.get_es_info(ene_method, thy_dct)
        ini_thy_info = filesys.inf.get_es_info(geo_model, thy_dct)
        pf_model = parser.model.pf_model_info(
            spc_model_dct[spc_model]['pf'])
        ts_dct[tsname] = parser.species.build_sing_chn_sadpt_dct(
            tsname, rxn, thy_info, ini_thy_info,
            run_inp_dct, spc_dct, cla_dct,
            direction='exo')
    spc_dct = parser.species.combine_sadpt_spc_dcts(
        ts_dct, spc_dct)

    # Build the MESS label idx dictionary for the PES
    label_dct = ktproutines.label.make_pes_label_dct(
        rxn_lst, pes_idx, spc_dct)

    # Set path where MESS files will be written and read
    mess_path = pfrunner.messrate_path(
        run_prefix, pes_formula, sub_pes_idx)

    # Try and read the MESS file from the filesystem first
    # _, _ = pfrunner.read_mess_file(mess_path)

    # Write the MESS file
    if write_messrate:  # and not mess_inp_str:
        print(('\n\n------------------------------------------------' +
               '--------------------------------------'))
        print('\nBuilding the MESS input file...')

        # Write the strings for the MESS input file
        globkey_str = ktproutines.rates.make_header_str(
            temps, pressures)

        # Write the energy transfer section strings for MESS file
        energy_trans_str = ktproutines.rates.make_etrans_str(
            spc_dct, rxn_lst, **etransfer)

        # Write the MESS strings for all the PES channels
        well_str, bi_str, ts_str, dats = ktproutines.rates.make_pes_mess_str(
            spc_dct, rxn_lst, pes_idx,
            run_prefix, save_prefix, label_dct,
            spc_model_dct, thy_dct)

        # Combine strings together
        mess_inp_str = ktproutines.rates.make_messrate_str(
            globkey_str, energy_trans_str,
            well_str, bi_str, ts_str)

        # Write the MESS file into the filesystem
        print(('\n++++++++++++++++++++++++++++++++++++++++++++++++' +
               '++++++++++++++++++++++++++++++++++++++'))
        print('\nWriting the MESS input file at {}'.format(mess_path))
        print(mess_inp_str)
        pfrunner.write_mess_file(mess_inp_str, dats, mess_path)

        # Write MESS file into job directory
        pfrunner.write_cwd_rate_file(mess_inp_str, pes_formula, sub_pes_idx)

    # Run mess to produce rate output
    if run_messrate:
        print(('\n\n------------------------------------------------' +
               '--------------------------------------'))
        print('\nRunning MESS for the input file at {}'.format(mess_path))
        pfrunner.run_rates(mess_path)

    # Fit rate output to modified Arrhenius forms, print in ChemKin format
    if run_fits:
        print(('\n\n------------------------------------------------' +
               '--------------------------------------'))
        print('\nFitting Rate Constants for PES to Functional Forms')
        ckin_str_lst = ktproutines.fit.fit_rates(
            temps, pressures, tunit, punit,
            pes_formula, label_dct,
            es_info, pf_model,
            mess_path, fit_method, pdep_fit,
            arrfit_thresh)
        writer.ckin.write_rates_file(ckin_rate_str_lst)
