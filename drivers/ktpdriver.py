""" driver for rate constant evaluations
"""

import autorun
import ratefit
from mechroutines.pf import ktp as ktproutines
from mechlib.amech_io import writer
from mechlib.amech_io import parser
from mechlib.amech_io import job_path
from mechlib.amech_io import output_path
from mechlib.amech_io import printer as ioprinter
from mechlib.reaction import split_unstable_rxn


def run(pes_rlst,
        ktp_tsk_lst,
        spc_dct, thy_dct,
        pes_model_dct, spc_model_dct,
        run_prefix, save_prefix):
    """ main driver for generation of full set of rate constants on a single PES
    """

    # --------------------------------------- #
    # LOOP OVER ALL OF THE SUBPES in PES_RLST #
    # --------------------------------------- #

    for pes_inf, rxn_lst in pes_rlst.items():

        # ---------------------------------------------- #
        # PREPARE INFORMATION TO PASS TO KTPDRIVER TASKS #
        # ---------------------------------------------- #

        # Set objects
        pes_formula, pes_idx, subpes_idx = pes_inf
        label_dct = None

        # Print PES Channels that are being run
        ioprinter.runlst(pes_inf, rxn_lst)

        # Set paths where files will be written and read
        mess_path = job_path(
            run_prefix, 'MESS', 'RATE', pes_formula, locs_idx=subpes_idx)

        # --------------------------------- #
        # RUN THE REQUESTED KTPDRIVER TASKS #
        # --------------------------------- #

        # Write the MESS file
        write_rate_tsk = parser.run.extract_task('write_mess', ktp_tsk_lst)
        if write_rate_tsk is not None:

            # Get all the info for the task
            tsk_key_dct = write_rate_tsk[-1]
            pes_model = tsk_key_dct['kin_model']
            spc_model = tsk_key_dct['spc_model']

            spc_dct, rxn_lst, label_dct = _process(
                pes_idx, rxn_lst, ktp_tsk_lst, spc_model_dct, spc_model,
                thy_dct, spc_dct, run_prefix, save_prefix)

            ioprinter.messpf('write_header')

            mess_inp_str, dats = ktproutines.rates.make_messrate_str(
                pes_idx, rxn_lst,
                pes_model, spc_model,
                spc_dct, thy_dct,
                pes_model_dct, spc_model_dct,
                label_dct,
                mess_path, run_prefix, save_prefix)

            autorun.write_input(
                mess_path, mess_inp_str,
                aux_dct=dats, input_name='mess.inp')

        # Run mess to produce rates (currently nothing from tsk lst keys used)
        run_rate_tsk = parser.run.extract_task('run_mess', ktp_tsk_lst)
        if run_rate_tsk is not None:

            ioprinter.obj('vspace')
            ioprinter.obj('line_dash')
            ioprinter.running('MESS for the input file', mess_path)
            autorun.run_script(autorun.SCRIPT_DCT['messrate'], mess_path)

        # Fit rate output to modified Arrhenius forms, print in ChemKin format
        run_fit_tsk = parser.run.extract_task('run_fits', ktp_tsk_lst)
        if run_fit_tsk is not None:

            # Get all the info for the task
            tsk_key_dct = run_fit_tsk[-1]
            pes_model = tsk_key_dct['kin_model']

            if label_dct is not None:
                spc_dct, rxn_lst, label_dct = _process(
                    pes_idx, rxn_lst, ktp_tsk_lst, spc_model_dct, spc_model,
                    thy_dct, spc_dct, run_prefix, save_prefix)

            ioprinter.obj('vspace')
            ioprinter.obj('line_dash')
            ioprinter.info_message(
                'Fitting Rate Constants for PES to Functional Forms',
                newline=1)

            ckin_str = ratefit.fit.fit_ktp_dct(
                mess_path=mess_path,
                inp_fit_method=pes_model_dct[pes_model]['fit_method'],
                pdep_dct=pes_model_dct[pes_model]['pdep_fit'],
                arrfit_dct=pes_model_dct[pes_model]['arrfit_fit'],
                chebfit_dct=pes_model_dct[pes_model]['chebfit_fit'],
                troefit_dct=pes_model_dct[pes_model]['troefit_fit'],
                label_dct=label_dct,
                fit_temps=pes_model_dct[pes_model]['rate_temps'],
                fit_pressures=pes_model_dct[pes_model]['pressures'],
                fit_tunit=pes_model_dct[pes_model]['tunit'],
                fit_punit=pes_model_dct[pes_model]['punit']
            )

            ckin_path = output_path('CKIN')
            writer.ckin.write_rxn_file(
                {pes_formula: ckin_str}, pes_formula, ckin_path)


# ------- #
# UTILITY #
# ------- #
def _process(pes_idx, rxn_lst, ktp_tsk_lst, spc_model_dct, spc_model,
             thy_dct, spc_dct, run_prefix, save_prefix):
    """ Build info needed for the task
    """

    spc_model_dct_i = spc_model_dct[spc_model]

    # Obtain all of the transitions states
    ioprinter.message(
        'Identifying reaction classes for transition states...')
    ts_dct = parser.spc.ts_dct_from_ktptsks(
        pes_idx, rxn_lst, ktp_tsk_lst, spc_model_dct, thy_dct,
        spc_dct, run_prefix, save_prefix)
    spc_dct = parser.spc.combine_sadpt_spc_dcts(
        ts_dct, spc_dct)

    # Set reaction list with unstable species broken apart
    ioprinter.message('Identifying stability of all species...', newline=1)
    chkd_rxn_lst = split_unstable_rxn(
        rxn_lst, spc_dct, spc_model_dct_i, thy_dct, save_prefix)

    # Build the MESS label idx dictionary for the PES
    print('chkd_rxn_lst', chkd_rxn_lst)
    label_dct = ktproutines.label.make_pes_label_dct(
        chkd_rxn_lst, pes_idx, spc_dct, spc_model_dct_i)

    return spc_dct, chkd_rxn_lst, label_dct
