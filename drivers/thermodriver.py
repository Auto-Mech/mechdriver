""" Driver for thermochemistry evaluations including
    heats-of-formation and NASA polynomials describing
    thermodynamic quantities: Enthalpy, Entropy, Gibbs
"""

import autorun
import mechanalyzer
import chemkin_io
import automol.inchi
from automol.inchi import formula_string as fstring
from mechanalyzer.inf import spc as sinfo
from mechroutines.pf import thermo as thmroutines
from mechroutines.pf import runner as pfrunner
from mechroutines.pf.models import ene
from mechlib import filesys
from mechlib.amech_io import writer
from mechlib.amech_io import parser
from mechlib.amech_io import printer as ioprinter
from mechlib.amech_io import thermo_paths
from mechlib.amech_io import job_path
from mechlib.amech_io import output_path
from mechlib.reaction import split_unstable_spc


def run(spc_rlst,
        therm_tsk_lst,
        pes_mod_dct, spc_mod_dct,
        spc_dct,
        run_prefix, save_prefix):
    """ main driver for thermo run
    """

    # Print Header fo
    ioprinter.info_message('Calculating Thermochem:')
    ioprinter.runlst(('SPC', 0, 0), spc_rlst)

    # ------------------------------------------------ #
    # PREPARE INFORMATION TO PASS TO THERMDRIVER TASKS #
    # ------------------------------------------------ #

    # Build a list of the species to calculate thermochem for loops below
    spc_mods = list(spc_mod_dct.keys())  # hack
    split_spc_lst = split_unstable_spc(
        spc_rlst, spc_dct, spc_mod_dct[spc_mods[0]], save_prefix)
    spc_queue = parser.rlst.spc_queue(
        'spc', tuple(split_spc_lst.values())[0])

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
                messpf_inp_str = thmroutines.qt.make_messpf_str(
                    pes_mod_dct[pes_mod]['therm_temps'],
                    spc_dct, spc_name,
                    pes_mod_dct[pes_mod], spc_mod_dct[spc_mod],
                    run_prefix, save_prefix)
                ioprinter.messpf('input_string')
                ioprinter.info_message(messpf_inp_str)
                autorun.write_input(
                    thm_paths[idx][spc_mod][0], messpf_inp_str,
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
                    pfrunner.mess.read_messpf(thm_paths[idx][spc_mod][0]))
            final_pf = pfrunner.mess.combine_pfs(_pfs, coeffs, operators)

            # need to clean thm path build
            tot_idx = len(spc_mods)
            spc_info = sinfo.from_dct(spc_dct[spc_name])
            spc_fml = automol.inchi.formula_string(spc_info[0])
            thm_prefix = [spc_fml, automol.inchi.inchi_key(spc_info[0])]
            thm_paths[idx]['final'] = (
                job_path(run_prefix, 'MESS', 'PF', thm_prefix, locs_idx=tot_idx),
                job_path(run_prefix, 'THERM', 'NASA', thm_prefix, locs_idx=tot_idx)
            )
            pfrunner.mess.write_mess_output(
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
            basis_dct, uniref_dct = thmroutines.basis.prepare_refs(
                ref_scheme, spc_dct, [[spc_name, None]],
                run_prefix, save_prefix)

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
            hf0k = thmroutines.heatform.calc_hform_0k(
                ene_spc, ene_basis, spc_basis, coeff_basis, ref_set=ref_enes)
            spc_dct[spc_name]['Hfs'] = [hf0k]

        # Write the NASA polynomials in CHEMKIN format
        ckin_nasa_str = ''
        ckin_path = output_path('CKIN')
        for idx, spc_name in enumerate(spc_queue):

            ioprinter.nasa('calculate', spc_name)

            # Write the header describing the models used in thermo calcs
            ckin_nasa_str += writer.ckin.model_header(spc_mods, spc_mod_dct)

            # Build and write the NASA polynomial in CHEMKIN-format string
            # Call dies if you haven't run "write mess" task
            ckin_nasa_str += thmroutines.nasapoly.build_polynomial(
                spc_name, spc_dct,
                thm_paths[idx]['final'][0], thm_paths[idx]['final'][1])
            ckin_nasa_str += '\n\n'

            print(ckin_nasa_str)
        nasa7_params_all =  chemkin_io.parser.thermo.create_spc_nasa7_dct(ckin_nasa_str)
        ioprinter.info_message('SPECIES\t\tH(0 K)[kcal/mol]\tH(298 K)[kcal/mol]\tS(298 K)[cal/mol K]\n')
        for spc_name in nasa7_params_all:
            nasa7_params =  nasa7_params_all[spc_name]
            h0 = spc_dct[spc_name]['Hfs'][0]
            h298 =  mechanalyzer.calculator.thermo.enthalpy(nasa7_params, 298.15) /1000.
            s298 =  mechanalyzer.calculator.thermo.entropy(nasa7_params, 298.15)
            ioprinter.info_message('{}\t{:3.2f}\t{:3.2f}\t{:3.2f}'.format(spc_name, h0, h298, s298))

        # Write all of the NASA polynomial strings
        writer.ckin.write_nasa_file(ckin_nasa_str, ckin_path)
