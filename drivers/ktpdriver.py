""" Driver for kinetic evaluations

    Main Loop of Driver:
        (1) PES

    Main Workflow:
        (1) Collate and process data from the SAVE filesystem
        (2) Format and write data into MESS input file
        (3) Run MESS file to obtain rate constants
        (5) Fit rate constants to functional forms
        (6) Write functional forms to mechanism file
"""

from mechroutines.ktp import tsk as ktp_tasks
from mechroutines.ktp import label as ktp_label
from mechlib.amech_io import parser
from mechlib.amech_io import rate_paths
from mechlib.amech_io import printer as ioprinter
from mechlib.reaction import split_unstable_pes
from ioformat import pathtools
from mess_io.reader import rates as mess_reader
from chemkin_io.writer import mechanism
from chemkin_io.writer import comments


def run(pes_rlst, pes_grp_dct,
        ktp_tsk_lst,
        spc_dct, glob_dct,
        pes_mod_dct, spc_mod_dct,
        run_prefix, save_prefix, mdriver_path):
    """ Executes all kinetics tasks.

        :param pes_rlst: species from PESs to run
        :type pes_rlst: tuple(dict[str: dict])
        :param spc_rlst: lst of species to run
        :type spc_rlst: tuple(dict[str: dict])
        :param es_tsk_lst: list of the electronic structure tasks
        :type es_tsk_lst: tuple(tuple(str, str, dict))
        :param spc_dct: species information
        :type spc_dct: dict[str:dict]
        :param glob_dct: global information for all species
        :type glob_dct: dict[str: dict]
        :param thy_dct: all of the theory information
        :type thy_dct: dict[str:dict]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
        :param mdriver_path: path where mechdriver is running
        :type mdriver_path: str
    """

    # ------------------------------------------------------------------ #
    # PREPARE GENERAL INFORMATION FOR ALL PES TO PASS TO KTPDRIVER TASKS #
    # ------------------------------------------------------------------ #

    write_rate_tsk = parser.run.extract_task('write_mess', ktp_tsk_lst)
    run_rate_tsk = parser.run.extract_task('run_mess', ktp_tsk_lst)
    run_fit_tsk = parser.run.extract_task('run_fits', ktp_tsk_lst)

    # Group the PESs into lists
    pes_grps_rlst = parser.rlst.pes_groups(pes_rlst, pes_grp_dct)

    # --------------------------------------- #
    # LOOP OVER ALL OF THE SUBPES in PES_RLST #
    # --------------------------------------- #

    print('pes groups')
    print(pes_grps_rlst)

    for (pes_grp_rlst, pes_param_dct) in pes_grps_rlst:

        print('WORKING ON PES GROUP NUM')
        print('\n\n')
        print(pes_grp_rlst)
        print(pes_param_dct)

        # Generate the paths needed for MESSRATE calculations
        rate_paths_dct = rate_paths(pes_grp_rlst, run_prefix)

        for pesgrp_num, (pes_inf, rxn_lst) in enumerate(pes_grp_rlst.items()):

            # ---------------------------------------------- #
            # PREPARE INFORMATION TO PASS TO KTPDRIVER TASKS #
            # ---------------------------------------------- #

            # Set objects
            _, pes_idx, _ = pes_inf
            label_dct = None

            # Print PES Channels that are being run
            ioprinter.runlst(pes_inf, rxn_lst)

            # --------------------------------- #
            # RUN THE REQUESTED KTPDRIVER TASKS #
            # --------------------------------- #

            # Write the MESS file
            if write_rate_tsk is not None:
                tsk_key_dct = write_rate_tsk[-1]
                spc_dct, rxn_lst, instab_chnls, label_dct = _process(
                    run_fit_tsk, ktp_tsk_lst, pes_idx, rxn_lst,
                    spc_mod_dct, spc_dct, glob_dct,
                    run_prefix, save_prefix)
                ktp_tasks.write_messrate_task(
                    pesgrp_num, pes_inf, rxn_lst,
                    tsk_key_dct, pes_param_dct,
                    spc_dct,
                    pes_mod_dct, spc_mod_dct,
                    instab_chnls, label_dct,
                    rate_paths_dct, run_prefix, save_prefix)

            # Run mess to produce rates (urrently nothing from tsk lst used)
            if run_rate_tsk is not None:
                ktp_tasks.run_messrate_task(rate_paths_dct, pes_inf)

            # Fit rates to functional forms; write parameters to ChemKin file
            if run_fit_tsk is not None:
                if label_dct is None:
                    spc_dct, rxn_lst, _, label_dct = _process(
                        run_fit_tsk, ktp_tsk_lst, pes_idx, rxn_lst,
                        spc_mod_dct, spc_dct, glob_dct,
                        run_prefix, save_prefix)
                ktp_tasks.run_fits_task(
                    pes_inf, rate_paths_dct, mdriver_path,
                    label_dct, pes_mod_dct, spc_mod_dct, tsk_key_dct)

                # Read MESS file and get rate constants
                mess_str = pathtools.read_file(mess_path, 'rate.out')
                rxn_ktp_dct = mess_reader.get_rxn_ktp_dct(
                    mess_str,
                    label_dct=label_dct,
                    filter_kts=True,
                    tmin=min(pes_mod_dct[pes_mod]['rate_temps']),
                    tmax=max(pes_mod_dct[pes_mod]['rate_temps']),
                    pmin=min(pes_mod_dct[pes_mod]['pressures']),
                    pmax=max(pes_mod_dct[pes_mod]['pressures'])
                )
             
                # Fit rates
                ratefit_dct = pes_mod_dct[pes_mod]['rate_fit']
                rxn_param_dct, rxn_err_dct = ratefit.fit.fit_rxn_ktp_dct(
                    rxn_ktp_dct,
                    ratefit_dct['fit_method'],
                    pdep_dct=ratefit_dct['pdep_fit'],
                    arrfit_dct=ratefit_dct['arrfit_fit'],
                    chebfit_dct=ratefit_dct['chebfit_fit'],
                    troefit_dct=ratefit_dct['troefit_fit'],
                )
            
                # Write the reactions block header, which contains model info
                rxn_block_cmt = writer.ckin.model_header((spc_mod,), spc_mod_dct)
    
                # Get the comments dct and write the Chemkin string
                rxn_cmts_dct = comments.get_rxn_cmts_dct(
                    rxn_err_dct=rxn_err_dct, rxn_block_cmt=rxn_block_cmt) 
                ckin_str = mechanism.write_chemkin_file(
                    rxn_param_dct=rxn_param_dct, rxn_cmts_dct=rxn_cmts_dct)
    
                # Write the file
                ckin_path = output_path('CKIN', prefix=mdriver_path)
                ckin_filename = pes_formula + '.ckin'
                pathtools.write_file(ckin_str, ckin_path, ckin_filename)


# ------- #
# UTILITY #
# ------- #
def _process(tsk, ktp_tsk_lst, pes_idx, rxn_lst,
             spc_mod_dct, spc_dct, glob_dct,
             run_prefix, save_prefix):
    """ Build info needed for the task
    """

    tsk_key_dct = tsk[-1]
    spc_mod = tsk_key_dct['spc_model']

    spc_mod_dct_i = spc_mod_dct[spc_mod]

    # Obtain all of the transitions states
    ioprinter.message(
        'Identifying reaction classes for transition states...')
    ts_dct = parser.spc.ts_dct_from_ktptsks(
        pes_idx, rxn_lst, ktp_tsk_lst, spc_mod_dct,
        spc_dct, run_prefix, save_prefix)
    spc_dct = parser.spc.combine_sadpt_spc_dcts(
        ts_dct, spc_dct, glob_dct)

    # Set reaction list with unstable species broken apart
    ioprinter.message('Identifying stability of all species...', newline=1)
    chkd_rxn_lst, instab_chnls = split_unstable_pes(
        rxn_lst, spc_dct, spc_mod_dct_i, save_prefix)

    # Build the MESS label idx dictionary for the PES
    label_dct = ktp_label.make_pes_label_dct(
        chkd_rxn_lst, pes_idx, spc_dct, spc_mod_dct_i)

    return spc_dct, chkd_rxn_lst, instab_chnls, label_dct
