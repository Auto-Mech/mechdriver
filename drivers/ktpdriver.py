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


def run(pes_rlst, pes_grp_dct,
        ktp_tsk_lst,
        spc_dct, glob_dct,
        thy_dct, pes_mod_dct, spc_mod_dct,
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

    for (pes_grp_rlst, pes_param_dct) in pes_grps_rlst:

        # print('WORKING ON PES GROUP NUM')
        # print(pes_grp_rlst)

        # Generate the paths needed for MESSRATE calculations
        rate_paths_dct = rate_paths(pes_grp_rlst, run_prefix)

        # Process info required ro run all of the PESs
        if write_rate_tsk is not None:
            proc_tsk = write_rate_tsk
        else:
            proc_tsk = run_fit_tsk
        spc_dct, all_rxn_lst, all_instab_chnls, label_dct = _process(
            proc_tsk, ktp_tsk_lst, pes_grp_rlst,
            spc_mod_dct, spc_dct, glob_dct,
            run_prefix, save_prefix)

        # ---------------------------------------- #
        # WRITE AND RUN TASK FOR EACH PES IN GROUP #
        # ---------------------------------------- #
        for pesgrp_num, (pes_inf, rxn_lst) in enumerate(pes_grp_rlst.items()):

            # Print PES Channels that are being run
            ioprinter.runlst(pes_inf, rxn_lst)

            # Write the MESS file
            if write_rate_tsk is not None:
                tsk_key_dct = write_rate_tsk[-1]
                ktp_tasks.write_messrate_task(
                    pesgrp_num, pes_inf, all_rxn_lst[pesgrp_num],
                    tsk_key_dct, pes_param_dct,
                    spc_dct,
                    thy_dct, pes_mod_dct, spc_mod_dct,
                    all_instab_chnls[pesgrp_num], label_dct,
                    rate_paths_dct, run_prefix, save_prefix)

            # Run mess to produce rates (urrently nothing from tsk lst used)
            if run_rate_tsk is not None:
                ktp_tasks.run_messrate_task(rate_paths_dct, pes_inf)

        # ---------------------------------------- #
        # FIT THE COMBINES RATES FOR ENTIRE GROUP  #
        # ---------------------------------------- #

        # Fit rates to functional forms; write parameters to ChemKin file
        if run_fit_tsk is not None:
            tsk_key_dct = run_fit_tsk[-1]
            ktp_tasks.run_fits_task(
                pes_grp_rlst, pes_param_dct, rate_paths_dct, mdriver_path,
                label_dct, pes_mod_dct, spc_mod_dct, thy_dct,
                tsk_key_dct)


# ------- #
# UTILITY #
# ------- #
def _process(tsk, ktp_tsk_lst, pes_grp_rlst,
             spc_mod_dct, spc_dct, glob_dct,
             run_prefix, save_prefix):
    """ Build info needed for the task
    """

    # Generic task/model info independent of PESs
    tsk_key_dct = tsk[-1]
    spc_mod = tsk_key_dct['spc_model']

    spc_mod_dct_i = spc_mod_dct[spc_mod]

    label_dct = {}
    all_chkd_rxn_lst, all_instab_chnls = (), ()
    for _, (pes_inf, rxn_lst) in enumerate(pes_grp_rlst.items()):

        _, pes_idx, _ = pes_inf

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

        all_chkd_rxn_lst += (chkd_rxn_lst,)
        all_instab_chnls += (instab_chnls,)

        # Build the MESS label idx dictionary for the PES
        label_dct.update(
            ktp_label.make_pes_label_dct(
                label_dct, chkd_rxn_lst, pes_idx,
                spc_dct, spc_mod_dct_i))

    return spc_dct, all_chkd_rxn_lst, all_instab_chnls, label_dct
