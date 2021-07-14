""" Driver for thermochemistry evaluations including
    heats-of-formation and NASA polynomials describing
    thermodynamic quantities: Enthalpy, Entropy, Gibbs

    Main Loop of Driver:
        (1) PES or SPC list

    Main Workflow:
        (1) Write MESS:
            (1) Collate and process data from the SAVE filesystem
            (2) Format and write data into MESSPF input file
        (2) Run MESS:
            (1) Run MESS file to obtain partition functions
        (3) Run fits:
            (1) Read all partition functions from MESS output
            (2) Run PACC+ThermP to fit thermo NASA polynomials
            (3) Write functional forms to mechanism file
"""

import autorun
import mechanalyzer
import chemkin_io
import automol.inchi
from automol.inchi import formula_string as fstring
from phydat import phycon
import thermfit
from mechanalyzer.inf import spc as sinfo
from mechroutines import thermo as thmroutines
from mechroutines.models import ene
from mechlib import filesys
from mechlib.amech_io import reader
from mechlib.amech_io import writer
from mechlib.amech_io import parser
from mechlib.amech_io import printer as ioprinter
from mechlib.amech_io import thermo_paths
from mechlib.amech_io import job_path
from mechlib.amech_io import output_path
from mechlib.reaction import split_unstable_full


def run(pes_rlst, spc_rlst,
        therm_tsk_lst,
        pes_mod_dct, spc_mod_dct,
        spc_dct,
        run_prefix, save_prefix, mdriver_path):
    """ Executes all thermochemistry tasks.

        :param pes_rlst: species from PESs to run
            [(PES formula, PES idx, SUP-PES idx)
            (CHANNEL idx, (REACS, PRODS))
        :type pes_rlst: tuple(dict[str: dict])
        :param spc_rlst: lst of species to run
        :type spc_rlst: tuple(dict[str: dict])
        :param es_tsk_lst: list of the electronic structure tasks
            tuple(tuple(obj, tsk, keyword_dict))
        :type es_tsk_lst: tuple(tuple(str, str, dict))
        :param spc_dct: species information
            dict[spc_name: spc_information]
        :type spc_dct: dict[str:dict]
        :param glob_dct: global information for all species
            dict[spc_name: spc_information]
        :type glob_dct: dict[str: dict]
        :param thy_dct: all of the theory information
            dict[thy name: inf]
        :type thy_dct: dict[str:dict]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
        :param mdriver_path: path where mechdriver is running
        :type mdriver_path: str
    """

    # Print Header
    ioprinter.info_message('Calculating Thermochem:')
    ioprinter.runlst(('SPC', 0, 0), spc_rlst)

    # ------------------------------------------------ #
    # PREPARE INFORMATION TO PASS TO THERMDRIVER TASKS #
    # ------------------------------------------------ #

    # Build a list of the species to calculate thermochem for loops below
    spc_mods = list(spc_mod_dct.keys())  # hack
    spc_mod_dct_i = spc_mod_dct[spc_mods[0]]
    split_rlst = split_unstable_full(
        pes_rlst, spc_rlst, spc_dct, spc_mod_dct_i, save_prefix)
    spc_queue = parser.rlst.spc_queue(
        tuple(split_rlst.values())[0], 'SPC')

    # Build the paths [(messpf, nasa)], models and levels for each spc
    thm_paths = thermo_paths(spc_dct, spc_queue, spc_mods, run_prefix)

    # ----------------------------------- #
    # RUN THE REQUESTED THERMDRIVER TASKS #
    # ----------------------------------- #

    # Write and Run MESSPF inputs to generate the partition functions
    write_messpf_tsk = parser.run.extract_task('write_mess', therm_tsk_lst)
    if write_messpf_tsk is not None:

        ioprinter.messpf('write_header')

        spc_mods, pes_mod = parser.models.extract_models(write_messpf_tsk)

        for idx, spc_name in enumerate(spc_queue):
            print('write test {}'.format(spc_name))
            for spc_mod in spc_mods:
                messpf_inp_str, dat_dct = thmroutines.qt.make_messpf_str(
                    pes_mod_dct[pes_mod]['therm_temps'],
                    spc_dct, spc_name,
                    pes_mod_dct[pes_mod], spc_mod_dct[spc_mod],
                    run_prefix, save_prefix)
                ioprinter.messpf('input_string')
                ioprinter.info_message(messpf_inp_str)
                autorun.write_input(
                    thm_paths[idx][spc_mod][0], messpf_inp_str,
                    aux_dct=dat_dct,
                    input_name='pf.inp')

    # Run the MESSPF files that have been written
    run_messpf_tsk = parser.run.extract_task('run_mess', therm_tsk_lst)
    if run_messpf_tsk is not None:

        spc_mod, pes_mod = parser.models.extract_models(run_messpf_tsk)
        spc_mods = parser.models.split_model(spc_mod[0])

        ioprinter.messpf('run_header')
        for idx, spc_name in enumerate(spc_queue):

            _spc_mods, coeffs, operators = spc_mods

            # Run MESSPF for all requested models, combine the PFS at the end
            ioprinter.message('Run MESSPF: {}'.format(spc_name), newline=1)
            _pfs = []
            for spc_mod in _spc_mods:
                autorun.run_script(
                   autorun.SCRIPT_DCT['messpf'],
                   thm_paths[idx][spc_mod][0])
                _pfs.append(
                    reader.mess.messpf(thm_paths[idx][spc_mod][0]))
            final_pf = thermfit.pf.combine(_pfs, coeffs, operators)

            # need to clean thm path build
            tdx = len(spc_mods)
            spc_info = sinfo.from_dct(spc_dct[spc_name])
            spc_fml = automol.inchi.formula_string(spc_info[0])
            thm_prefix = [spc_fml, automol.inchi.inchi_key(spc_info[0])]
            thm_paths[idx]['final'] = (
                job_path(run_prefix, 'MESS', 'PF', thm_prefix, locs_idx=tdx),
                job_path(run_prefix, 'THERM', 'NASA', thm_prefix, locs_idx=tdx)
            )
            writer.mess.output(
                fstring(spc_dct[spc_name]['inchi']),
                final_pf, thm_paths[idx]['final'][0],
                filename='pf.dat')

    # Use MESS partition functions to compute thermo quantities
    run_fit_tsk = parser.run.extract_task('run_fits', therm_tsk_lst)
    if run_fit_tsk is not None:

        spc_mods, pes_mod = parser.models.extract_models(run_fit_tsk)
        pes_mod_dct_i = pes_mod_dct[pes_mod]

        ioprinter.nasa('header')
        chn_basis_ene_dct = {}
        for idx, spc_name in enumerate(spc_queue):

            # Take species model and add it to the chn_basis_ene dct
            spc_mod = spc_mods[0]
            spc_mod_dct_i = spc_mod_dct[spc_mod]
            if spc_mod not in chn_basis_ene_dct:
                chn_basis_ene_dct[spc_mod] = {}

            # Get the reference scheme and energies (ref in different place)
            ref_scheme = pes_mod_dct_i['therm_fit']['ref_scheme']
            ref_enes = pes_mod_dct_i['therm_fit']['ref_enes']

            # Determine info about the basis species used in thermochem calcs
            basis_dct, uniref_dct = thermfit.prepare_refs(
                ref_scheme, spc_dct, (spc_name,))

            # Get the basis info for the spc of interest
            spc_basis, coeff_basis = basis_dct[spc_name]

            # Get the energies for the spc and its basis
            ene_basis = []
            energy_missing = False
            for spc_basis_i in spc_basis:
                if spc_basis_i in chn_basis_ene_dct[spc_mod]:
                    ioprinter.message(
                        'Energy already found for basis species: '
                        + spc_basis_i)
                    ene_basis.append(chn_basis_ene_dct[spc_mod][spc_basis_i])
                else:
                    ioprinter.message(
                        'Energy will be determined for basis species: '
                        + spc_basis_i)
                    energy_missing = True
            if not energy_missing:
                pf_filesystems = filesys.models.pf_filesys(
                    spc_dct[spc_name], spc_mod_dct_i,
                    run_prefix, save_prefix, saddle=False)
                ene_spc = ene.read_energy(
                    spc_dct[spc_name], pf_filesystems, spc_mod_dct_i,
                    run_prefix, read_ene=True, read_zpe=True, saddle=False)
            else:
                ene_spc, ene_basis = thmroutines.basis.basis_energy(
                    spc_name, spc_basis, uniref_dct, spc_dct,
                    spc_mod_dct_i,
                    run_prefix, save_prefix)
                for spc_basis_i, ene_basis_i in zip(spc_basis, ene_basis):
                    chn_basis_ene_dct[spc_mod][spc_basis_i] = ene_basis_i

            # Calculate and store the 0 K Enthalpy
            hf0k = thermfit.heatform.calc_hform_0k(
                ene_spc, ene_basis, spc_basis, coeff_basis, ref_set=ref_enes)
            spc_dct[spc_name]['Hfs'] = [hf0k]

        # Write the NASA polynomials in CHEMKIN format
        ckin_nasa_str = ''
        ckin_path = output_path('CKIN', prefix=mdriver_path)
        for idx, spc_name in enumerate(spc_queue):

            ioprinter.nasa('calculate', spc_name)

            # Write the header describing the models used in thermo calcs
            ckin_nasa_str += writer.ckin.model_header(
                spc_mods, spc_mod_dct, refscheme=ref_scheme)

            # Build and write the NASA polynomial in CHEMKIN-format string
            # Call dies if you haven't run "write mess" task
            print(thm_paths[idx])
            ckin_nasa_str += thmroutines.nasapoly.build_polynomial(
                spc_name, spc_dct,
                thm_paths[idx]['final'][0], thm_paths[idx]['final'][1])
            ckin_nasa_str += '\n\n'
        print('CKIN NASA STR\n')
        print(ckin_nasa_str)

        nasa7_params_all = chemkin_io.parser.thermo.create_spc_nasa7_dct(ckin_nasa_str)
        # print('ckin_nasa_str test', ckin_nasa_str)
        templist = (298.15, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500)
        for spc_name in nasa7_params_all:
            nasa7_params = nasa7_params_all[spc_name]
            whitespace = 18-len(spc_name)
            whitespace = whitespace*' '
            hf0 = spc_dct[spc_name]['Hfs'][0] * phycon.EH2KCAL
            hf298 = mechanalyzer.calculator.thermo.enthalpy(nasa7_params, 298.15) /1000.
            ioprinter.info_message('SPECIES            H0f(0 K)  H0f(298 K) in kcal/mol:')
            ioprinter.info_message('{}{}{:>9.2f}{:>9.2f}'
                                   .format(spc_name, whitespace, hf0, hf298))
            ioprinter.info_message('\n T (K)   H - H(T)    S(T)      Cp(T) ')
            ioprinter.info_message('Kelvin  kcal/mol cal/(mol K) cal/(mol K)')
            hincref = hf298
            # hincref = hf298 - (mechanalyzer.calculator.thermo.enthalpy(nasa7_params, 10) /1000.)
            for temp in templist:
                hinct = mechanalyzer.calculator.thermo.enthalpy(nasa7_params, temp) /1000. - hincref
                entt = mechanalyzer.calculator.thermo.entropy(nasa7_params, temp)
                cpt = mechanalyzer.calculator.thermo.heat_capacity(nasa7_params, temp)
                ioprinter.info_message('{:>7.2f}{:>9.2f}{:>9.2f}{:>9.2f}'
                                       .format(temp, hinct, entt, cpt))

        # Write all of the NASA polynomial strings
        writer.ckin.write_nasa_file(ckin_nasa_str, ckin_path)
