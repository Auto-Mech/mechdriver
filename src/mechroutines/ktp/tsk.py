""" Tasks for kTPDriver
"""

import os
import ioformat
import chemkin_io
import autorun
import ratefit
from mechlib import filesys
from mechlib.amech_io import writer
from mechlib.amech_io import output_path
from mechlib.amech_io import printer as ioprinter
from mechroutines.models.typ import is_abstraction_pes
from mechroutines.ktp.rates import make_full_str
from mechroutines.ktp.rates import make_global_etrans_str
from mechroutines.ktp.rates import make_pes_mess_str
from mechroutines.ktp._multipes import obtain_multipes_rxn_ktp_dct


def write_messrate_task(pesgrp_num, pes_inf, rxn_lst,
                        tsk_key_dct, pes_param_dct,
                        spc_dct,
                        thy_dct, pes_model_dct, spc_model_dct,
                        unstab_chnls, label_dct,
                        rate_paths_dct, run_prefix, save_prefix,
                        nprocs=1):
    """ Reads and processes all information in the save filesys for
        all species on the PES that are required for MESS rate calculations,
        as specified by the model dictionaries built from user input.

        :param pes_idx:
        :type pes_idx: int
        :param rxn_lst:
        :type rxn_lst:
        :param pes_model: model for PES conditions for rates from user input
        :type pes_model: str
        :param spc_model: model for partition fxns for rates from user input
        :type spc_model: str
        :param mess_path: path to write mess file (change since pfx given?)
    """

    _, pes_idx, _ = pes_inf

    pes_mod = tsk_key_dct['kin_model']
    spc_mod = tsk_key_dct['spc_model']

    pes_model_dct_i = pes_model_dct[pes_mod]
    spc_model_dct_i = spc_model_dct[spc_mod]

    # Write the MESS strings for all the PES channels
    rxn_chan_str, dats, hot_enes_dct = make_pes_mess_str(
        spc_dct, rxn_lst, pes_idx, pesgrp_num, unstab_chnls,
        run_prefix, save_prefix, label_dct,
        tsk_key_dct, pes_param_dct,
        thy_dct, pes_model_dct_i, spc_model_dct_i, spc_mod,
        nprocs=nprocs)

    # Write the energy transfer section strings for MESS file
    energy_trans_str = make_global_etrans_str(
        rxn_lst, spc_dct, pes_model_dct_i['energy_transfer'])

    # Create full string by writing the appropriate header, accounting for
    # (1) MESS Version and (2) Use of Well-Extension
    # And include the global_etrans and reaction channel strings
    pes_param_dct = make_full_str(
        energy_trans_str, rxn_chan_str, dats,
        pesgrp_num, pes_param_dct, hot_enes_dct,
        rate_paths_dct, pes_inf,
        pes_model_dct_i,
        spc_dct, rxn_lst, pes_idx, tsk_key_dct)

    return pes_param_dct


def run_messrate_task(pes_inf, rxn_lst, tsk_key_dct, spc_dct, rate_paths_dct):
    """ Run the MESSRATE input file.

        First tries to run a well-lumped file, then tries to
        run the base file if it exists.

        Need an overwrite task
    """

    _, pes_idx, _ = pes_inf

    # Get the path to the MESSRATE file to run
    # (1) Vers1-Base, (2) Vers1-WellLump, (3) Vers2-Base
    path_dct = rate_paths_dct[pes_inf]
    mess_version = tsk_key_dct['mess_version']
    if (
        mess_version == 'v1' and
        tsk_key_dct['well_lumping'] and
        not is_abstraction_pes(spc_dct, rxn_lst, pes_idx)
    ):
        path = path_dct[f'wext-{mess_version}']
        typ = 'wext'
    else:
        path = path_dct[f'base-{mess_version}']
        typ = 'base'

    mess_inp = os.path.join(path, 'mess.inp')
    # mess_out = os.path.join(path, 'mess.out')
    # if not os.path.exists(mess_out): (rerun option needed?)
    if os.path.exists(mess_inp):
        ioprinter.obj('vspace')
        ioprinter.obj('line_dash')
        if typ == 'base':
            ioprinter.running(
                f'MESS base input with version {mess_version} '
                f'at {path}')
        else:
            ioprinter.running(
                f'MESS well-lumped input with version {mess_version} '
                f'at {path}')
        autorun.run_script(
            autorun.SCRIPT_DCT[f'messrate-{mess_version}'], path
        )
    else:
        if typ == 'base':
            ioprinter.warning_message(
                f'No MESS base input for version {mess_version} '
                f'found at {path}')
        else:
            ioprinter.warning_message(
                f'No MESS well-lumped input for version {mess_version} '
                f'found at {path}')


def run_fits_task(pes_grp_rlst, pes_param_dct, rate_paths_dct, mdriver_path,
                  pes_mod_dct, spc_mod_dct, thy_dct,
                  tsk_key_dct, spc_dct):
    """ Run the fits and potentially

        assume that the rate_paths_dct will come in with all PESs in group
    """

    # Combine all PESs into a string for writing the CKIN file
    pes_strs = ()
    for pes_inf in pes_grp_rlst.keys():
        _inf = (pes_inf[0], str(pes_inf[1]+1), str(pes_inf[2]+1))
        pes_strs += ('_'.join(_inf),)
    tot_fml = '-'.join(pes_strs)

    # Get the model and sort info from tsk key dct
    pes_mod = tsk_key_dct['kin_model']
    spc_mod = tsk_key_dct['spc_model']
    sort_info_lst = filesys.mincnf.sort_info_lst(tsk_key_dct['sort'], thy_dct)

    ioprinter.obj('vspace')
    ioprinter.obj('line_dash')

    # Obtain the rate constants from the MESS files
    ioprinter.info_message(
        'Reading Rate Constants from MESS outputs', newline=1)
    rxn_ktp_dct = obtain_multipes_rxn_ktp_dct(
        pes_grp_rlst, rate_paths_dct, pes_param_dct,
        pes_mod_dct, pes_mod,
        tsk_key_dct, spc_dct)

    # Fit the rate constants
    ioprinter.info_message(
        'Fitting Rate Constants for PES to Functional Forms', newline=1)
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
    rxn_block_cmt = writer.ckin.model_header(
        (spc_mod,), spc_mod_dct,
        sort_info_lst=sort_info_lst,
        refscheme=pes_mod_dct[pes_mod]['therm_fit']['ref_scheme'])

    # Get the comments dct and write the Chemkin string
    rxn_cmts_dct = chemkin_io.writer.comments.get_rxn_cmts_dct(
        rxn_err_dct=rxn_err_dct, rxn_block_cmt=rxn_block_cmt)
    ckin_str = chemkin_io.writer.mechanism.write_chemkin_file(
        rxn_param_dct=rxn_param_dct, rxn_cmts_dct=rxn_cmts_dct)

    # Write the file
    ckin_path = output_path('CKIN', prefix=mdriver_path)
    ckin_filename = f'{tot_fml}.ckin'
    ioformat.pathtools.write_file(ckin_str, ckin_path, ckin_filename)
