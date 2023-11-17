""" Runs automech instancs for tests
"""

import os
import shutil
import chemkin_io
from ioformat import pathtools
import mechanalyzer
import ratefit
from mechlib.filesys import prefix_fs
from mechlib.amech_io import parser as ioparser
from mechlib.amech_io import printer as ioprinter
from drivers import esdriver, thermodriver, ktpdriver, transdriver, procdriver


# Helper functions to run a single instance of MechDriver
def run_mechdriver(run_template,
                   tmp_dir,
                   tmp_inp_dir, cwd_inp_dir,
                   tmp_run_dir, tmp_save_dir):
    """ test automech.py
    """
    # Copy input to tmp directory and replace it
    if not os.path.exists(tmp_inp_dir):
        shutil.copytree(cwd_inp_dir, tmp_inp_dir)

    print('TMP DIR PRINT')
    print('inp', tmp_inp_dir)
    print('run', tmp_run_dir)
    print('save', tmp_save_dir)
    print('run template', run_template)

    _fill_template_and_write_file(
        run_template, 'run.dat',
        tmp_inp_dir, tmp_run_dir, tmp_save_dir)

    _mechdriver_main(tmp_dir)


def _fill_template_and_write_file(templatefile, inpfile,
                                  tmp_inp_dir, tmp_run_dir, tmp_save_dir):
    """ Read the run.dat and replace run_prefix and save_prefix
    """

    # Set names of template and input for the calculation
    inp_file = os.path.join(tmp_inp_dir, templatefile)
    new_inp_file = os.path.join(tmp_inp_dir, inpfile)

    # Read template and fill with the run and save prefix
    with open(inp_file, 'r') as fobj:
        inp_str = fobj.read()
    new_inp_str = inp_str.format(tmp_run_dir, tmp_save_dir)

    # Write the run.dat and models.dat for the calculation
    with open(new_inp_file, 'w') as fobj:
        fobj.write(new_inp_str)


def _mechdriver_main(tmp_dir):
    """ Copy of MechDriver bin
    """

    # print header message and host name (probably combine into one function)
    ioprinter.program_header('amech')
    ioprinter.random_cute_animal()
    ioprinter.host_name()

    # parse all of the input
    ioprinter.program_header('inp')

    inp_strs = ioparser.read_amech_input(tmp_dir)

    thy_dct = ioparser.thy.theory_dictionary(inp_strs['thy'])
    kmod_dct, smod_dct = ioparser.models.models_dictionary(
        inp_strs['mod'], thy_dct)
    inp_key_dct = ioparser.run.input_dictionary(inp_strs['run'])
    pes_idx_dct, spc_idx_dct = ioparser.run.chem_idxs(inp_strs['run'])
    tsk_lst_dct = ioparser.run.tasks(
        inp_strs['run'],  thy_dct)
        #inp_strs['run'], inp_strs['mech'], thy_dct)
    spc_dct, glob_dct = ioparser.spc.species_dictionary(
        inp_strs['spc'], inp_strs['dat'], inp_strs['geo'], {}, inp_key_dct, 'csv')
    pes_dct = ioparser.mech.pes_dictionary(
        inp_strs['mech'], 'chemkin', spc_dct)

    pes_rlst, spc_rlst = ioparser.rlst.run_lst(
        pes_dct, spc_dct, pes_idx_dct, spc_idx_dct)

    # build the run-save filesystem directories
    prefix_fs(inp_key_dct['run_prefix'], inp_key_dct['save_prefix'])

    # run drivers requested by user
    es_tsks = tsk_lst_dct.get('es')
    if es_tsks is not None:
        ioprinter.program_header('es')
        esdriver.run(
            pes_rlst, spc_rlst,
            es_tsks,
            spc_dct, glob_dct, thy_dct,
            inp_key_dct['run_prefix'], inp_key_dct['save_prefix']
        )
        ioprinter.program_exit('es')

    therm_tsks = tsk_lst_dct.get('thermo')
    if therm_tsks is not None:
        ioprinter.program_header('thermo')
        thermodriver.run(
            pes_rlst, spc_rlst,
            therm_tsks,
            kmod_dct, smod_dct,
            spc_dct, thy_dct,
            inp_key_dct['run_prefix'], inp_key_dct['save_prefix'], tmp_dir
        )
        ioprinter.program_exit('thermo')

    trans_tsks = tsk_lst_dct.get('trans')
    if trans_tsks is not None:
        ioprinter.program_header('trans')
        if pes_dct:
            transdriver.run(
                pes_rlst, spc_rlst,
                trans_tsks,
                smod_dct,
                spc_dct, thy_dct,
                inp_key_dct['run_prefix'], inp_key_dct['save_prefix']
            )
        ioprinter.program_exit('trans')

    ktp_tsks = tsk_lst_dct.get('ktp')
    if ktp_tsks is not None:
        ioprinter.program_header('ktp')
        ktpdriver.run(
            pes_rlst,
            ktp_tsks,
            spc_dct, glob_dct,
            kmod_dct, smod_dct,
            inp_key_dct['run_prefix'], inp_key_dct['save_prefix'], tmp_dir
        )
        ioprinter.program_exit('ktp')

    proc_tsks = tsk_lst_dct.get('proc')
    if proc_tsks is not None:
        ioprinter.program_header('proc')
        procdriver.run(
            pes_rlst, spc_rlst,
            proc_tsks,
            spc_dct, thy_dct,
            kmod_dct, smod_dct,
            inp_key_dct['run_prefix'], inp_key_dct['save_prefix'], tmp_dir
        )
        ioprinter.program_exit('proc')

    # exit program
    ioprinter.obj('vspace')
    ioprinter.program_exit('amech')


# Helper functions to ensure consistent numerical results
def chk_therm(therm_dat_file, therm_calc_file,
              therm_dat_path, therm_calc_path,
              tmp_inp_dir,
              temps):
    """ Read Check the
    """

    # Read the data in the thermo and rate CKIN files
    ckin_path = os.path.join(therm_calc_path, 'CKIN')

    therm_calc_str = pathtools.read_file(ckin_path, therm_calc_file)
    therm_dat_str = pathtools.read_file(therm_dat_path, therm_dat_file)

    nasa7_calc = chemkin_io.parser.thermo.create_spc_nasa7_dct(therm_calc_str)
    nasa7_dat = chemkin_io.parser.thermo.create_spc_nasa7_dct(therm_dat_str)

    thm_calc = mechanalyzer.calculator.thermo.create_spc_thermo_dct(
        nasa7_calc, temps)
    thm_dat = mechanalyzer.calculator.thermo.create_spc_thermo_dct(
        nasa7_dat, temps)

    spc_str = pathtools.read_file(tmp_inp_dir, 'species.csv')
    spc_ident_dct = mechanalyzer.parser.spc.build_spc_dct(spc_str, 'csv')

    thm_dct = mechanalyzer.calculator.compare.get_aligned_spc_thermo_dct(
        [thm_dat, thm_calc], [spc_ident_dct, spc_ident_dct])

    for thm_data in thm_dct.values():
        # 0 is data/therm and 1 is calc'd therm
        assert _assess(thm_data[0][0], thm_data[1][0], thresh=3.0)
        assert _assess(thm_data[0][1], thm_data[1][1], thresh=3.0)
        assert _assess(thm_data[0][2], thm_data[1][2], thresh=3.0)
        assert _assess(thm_data[0][3], thm_data[1][3], thresh=3.0)


def chk_rates(rates_dat_file, rates_calc_file,
              rates_dat_path,
              tmp_dir,  # tmp_inp_dir,
              pressures, temps):
    """ Read Check the
    """

    # Read the data in the thermo and rate CKIN files
    ckin_path = os.path.join(tmp_dir, 'CKIN')

    rates_calc_str = pathtools.read_file(ckin_path, rates_calc_file)
    rates_dat_str = pathtools.read_file(rates_dat_path, rates_dat_file)

    rxn_str_calc = chemkin_io.parser.mechanism.reaction_block(rates_calc_str)
    rxn_str_dat = chemkin_io.parser.mechanism.reaction_block(rates_dat_str)
    units = chemkin_io.parser.mechanism.reaction_units(rates_dat_str)

    par_dct_calc = chemkin_io.parser.reaction.param_dct(rxn_str_calc, *units)
    par_dct_dat = chemkin_io.parser.reaction.param_dct(rxn_str_dat, *units)

    rxn_ktp_dct_calc = mechanalyzer.calculator.rates.eval_rxn_param_dct(
        par_dct_calc, pressures, temps)
    rxn_ktp_dct_dat = mechanalyzer.calculator.rates.eval_rxn_param_dct(
        par_dct_dat, pressures, temps)

    # spc_str = pathtools.read_file([], 'species.csv', path=TMP_INP_DIR)
    # spc_ident_dct = mechanalyzer.parser.spc.build_spc_dct(spc_str, 'csv')

    # ktp_dct = mechanalyzer.calculator.compare.get_aligned_rxn_ktp_dct(
    #     [ktp_dct_dat, ktp_dct_calc], [],
    #     [spc_ident_dct, spc_ident_dct],
    #     (500., 1000., 1500., 2000.))
    # print(ktp_dct)

    for rxn, ktp_dct_calc in rxn_ktp_dct_calc.items():
        for pressure in ktp_dct_calc:
            calc = ratefit.ktpdct.read(
                rxn_ktp_dct_calc, rxn, pressure, 'rates')
            dat = ratefit.ktpdct.read(
                rxn_ktp_dct_dat, rxn, pressure, 'rates')
            assert _assess(dat, calc, thresh=3.0)


def _assess(dat1, dat2, thresh):
    """ Assess if two sets of values are within some numerical threshold
    """

    cond = True
    for val1, val2 in zip(dat1, dat2):
        cond = bool((abs(val1 - val2) / val1) * 100.0 < thresh)

    return cond
