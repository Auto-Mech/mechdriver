""" driver for rate constant evaluations
"""

import os
import autorun
from mechanalyzer.inf import thy as tinfo
from mechroutines.pf import ktp as ktproutines
from mechroutines.pf import runner as pfrunner
from mechlib import filesys
from mechlib.amech_io import writer
from mechlib.amech_io import parser
from mechlib.amech_io import printer as ioprinter
from mechlib.structure import instab


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
        pes_model_dct[pes_model]['dbl_arrfit_check']
    )

    # Obtain all of the transitions states
    ioprinter.message('Identifitying reaction classes for transition states...')
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
        method_dct = thy_dct.get(ene_method)
        ini_method_dct = thy_dct.get(geo_model)
        thy_info = tinfo.from_dct(method_dct)
        ini_thy_info = tinfo.from_dct(ini_method_dct)
        pf_model = parser.model.pf_model_info(
            spc_model_dct[spc_model]['pf'])
        ts_dct[tsname] = parser.species.build_sing_chn_sadpt_dct(
            tsname, rxn, thy_info, ini_thy_info,
            run_inp_dct, spc_dct, cla_dct, run_prefix, save_prefix,
            direction='forw')
    spc_dct = parser.species.combine_sadpt_spc_dcts(
        ts_dct, spc_dct)

    # Set reaction list with unstable species broken apart
    ioprinter.message('Identifying stability of all species...', newline=1)
    # rxn_lst = instab.break_all_unstable(
    #     rxn_lst, spc_dct, spc_model_dct, thy_dct, save_prefix)
    # Build the MESS label idx dictionary for the PES
    label_dct = ktproutines.label.make_pes_label_dct(
        rxn_lst, pes_idx, spc_dct, spc_model_dct)

    # Set paths where files will be written and read
    mess_path = pfrunner.messrate_path(
        run_prefix, pes_formula, sub_pes_idx)
    starting_path = os.getcwd()
    ckin_path = os.path.join(starting_path, 'ckin')

    # Write the MESS file
    if write_messrate:  # and not mess_inp_str:
        
        ioprinter.messpf('write_header')

        # Write the strings for the MESS input file
        globkey_str = ktproutines.rates.make_header_str(
            temps, pressures)

        # Write the energy transfer section strings for MESS file
        energy_trans_str = ktproutines.rates.make_global_etrans_str(
            rxn_lst, spc_dct, etransfer)

        # Write the MESS strings for all the PES channels
        chan_str, dats, p_enes, cnlst = ktproutines.rates.make_pes_mess_str(
            spc_dct, rxn_lst, pes_idx,
            run_prefix, save_prefix, label_dct,
            spc_model_dct, thy_dct)

        # Combine strings together
        mess_inp_str = ktproutines.rates.make_messrate_str(
            globkey_str, energy_trans_str, chan_str)

        # Write the MESS file into the filesystem
        ioprinter.obj('line_plus')
        ioprinter.writing('MESS input file', mess_path)
        ioprinter.debug_message(mess_inp_str)

        autorun.write_input(
            mess_path, mess_inp_str,
            aux_dct=dats, input_name='mess.inp')

        # Write MESS file into job directory
        pfrunner.write_cwd_rate_file(mess_inp_str, pes_formula, sub_pes_idx)

        # Create a plot of the PES energies (not working correctly)
        # ktproutines.plot_from_dct(p_enes, cnlst, pes_formula)

    # Run mess to produce rate output
    if run_messrate:
        ioprinter.obj('vspace')
        ioprinter.obj('line_dash')
        ioprinter.running('MESS for the input file', mess_path)
        autorun.run_script(autorun.SCRIPT_DCT['messrate'], mess_path)

    # Fit rate output to modified Arrhenius forms, print in ChemKin format
    if run_fits:
        ioprinter.obj('vspace')
        ioprinter.obj('line_dash')
        ioprinter.info_message('Fitting Rate Constants for PES to Functional Forms', newline=1)
        ckin_str_dct = ktproutines.fit.fit_rates(
            temps, pressures, tunit, punit,
            pes_formula, label_dct,
            es_info, pf_model,
            mess_path, fit_method, pdep_fit,
            arrfit_thresh)
        writer.ckin.write_rxn_file(ckin_str_dct, pes_formula, ckin_path)
