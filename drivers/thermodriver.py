""" Driver for thermochemistry evaluations including
    heats-of-formation and NASA polynomials describing
    thermodynamic quantities: Enthalpy, Entropy, Gibbs
"""

import autorun
import automol.inchi
from mechroutines.pf import thermo as thmroutines
from mechroutines.pf import runner as pfrunner
from mechroutines.pf.models import ene
from mechlib.amech_io import writer
from mechlib.amech_io import parser
from mechlib.amech_io import printer as ioprinter
from mechlib.amech_io import thermo_paths
from mechlib.amech_io import job_path
from mechlib.amech_io import output_path
from mechanalyzer.inf import spc as sinfo
from mechlib.reaction import split_unstable_spc
from mechlib import filesys
from automol.inchi import formula_string as fstring


def run(spc_rlst,
        therm_tsk_lst,
        pes_model_dct, spc_model_dct,
        spc_dct, thy_dct,
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
    split_spc_lst = split_unstable_spc(
        spc_rlst, spc_dct, spc_model_dct, thy_dct, save_prefix)
    spc_queue = parser.species.build_queue(split_spc_lst)

    # Build the paths [(messpf, nasa)], models and levels for each spc
    thm_paths = thermo_paths(spc_dct, spc_queue, run_prefix)

    # Get models from the task
    pf_levels, pf_models = parser.models.mult_models(
        mods, spc_model_dct, thy_dct)

    # ----------------------------------- #
    # RUN THE REQUESTED THERMDRIVER TASKS #
    # ----------------------------------- #

    # for spc in spc_rlst.values():

    # Write and Run MESSPF inputs to generate the partition functions
    write_messpf_tsk = parser.run.extract_tsk('write_messpf', therm_tsk_lst)
    if write_messpf_tsk is not None:

        spc_models, pes_model = parser.run.extract_models(write_messpf_tsk)

        ioprinter.messpf('write_header')
        for idx, spc_name in enumerate(spc_queue):
            for spc_model in spc_models:
                messpf_inp_str = thmroutines.qt.make_messpf_str(
                    pes_model_dct[pes_model]['therm_temps'],
                    spc_dct, spc_name,
                    pf_models[spc_model], pf_levels[spc_model],
                    run_prefix, save_prefix)
                ioprinter.messpf('input_string')
                autorun.write_input(
                    thm_paths[idx][spc_model][0], messpf_inp_str,
                    input_name='pf.inp')

    # Run the MESSPF files that have been written
    run_messpf_tsk = parser.tsks.extract_tsk('run_messpf', therm_tsk_lst)
    if run_messpf_tsk is not None:

        spc_model, pes_model = parser.run.extract_models(run_messpf_tsk)
        spc_models = parser.models.split_model(spc_model)

        ioprinter.messpf('run_header')
        for idx, spc_name in enumerate(spc_queue):

            _spc_models, coeffs, operators = spc_models

            # Run MESSPF for all requested models, combine the PFS at the end
            ioprinter.message('{}'.format(spc_name), newline=1)
            _pfs = []
            for spc_mod in _spc_models:
                autorun.run_script(
                   autorun.SCRIPT_DCT['messpf'],
                   thm_paths[idx][spc_mod][0])
                _pfs.append(pfrunner.mess.read_messpf(
                    thm_paths[idx][spc_mod][0]))
            final_pf = pfrunner.mess.combine_pfs(_pfs, coeffs, operators)

            # need to clean thm path build
            totidx = len(spc_models)
            spc_info = sinfo.from_dct(spc_dct[spc_name])
            spc_fml = automol.inchi.formula_string(spc_info[0])
            thm_paths[idx]['final'] = (
                job_path(run_prefix, 'MESS', 'PF', spc_fml, locs_idx=totidx),
                job_path(run_prefix, 'THERM', 'NASA', spc_fml, locs_idx=totidx)
            )
            pfrunner.mess.write_mess_output(
                fstring(spc_dct[spc_name]['inchi']),
                final_pf, thm_paths[idx]['final'][0],
                filename='pf.dat')

    # Use MESS partition functions to compute thermo quantities
    run_fit_tsk = parser.tsks.extract_tsk('run_fit', therm_tsk_lst)
    if run_fit_tsk is not None:

        spc_models, pes_model = parser.run.extract_models(run_messpf_tsk)

        ioprinter.nasa('header')
        chn_basis_ene_dct = {}
        for idx, spc_name in enumerate(spc_queue):

            # Take species model and add it to the chn_basis_ene dct
            spc_model = spc_models[0]
            if spc_model not in chn_basis_ene_dct:
                chn_basis_ene_dct[spc_model] = {}

            # Get the reference scheme and energies (ref in different place)
            ref_scheme = spc_model_dct[spc_model]['therm_fit']['ref_scheme']
            ref_enes = spc_model_dct[spc_model]['therm_fit']['ref_enes']

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
                if spc_basis_i in chn_basis_ene_dct[spc_model]:
                    ioprinter.message(
                        'Energy already found for basis species: '
                        + spc_basis_i)
                    ene_basis.append(chn_basis_ene_dct[spc_model][spc_basis_i])
                else:
                    ioprinter.message(
                        'Energy will be determined for basis species: '
                        + spc_basis_i)
                    energy_missing = True
            if not energy_missing:
                pf_filesystems = filesys.models.pf_filesys(
                    spc_dct[spc_name], pf_levels[spc_model],
                    run_prefix, save_prefix, saddle=False)
                ene_spc = ene.read_energy(
                    spc_dct[spc_name], pf_filesystems, pf_models[spc_model],
                    pf_levels[spc_model],
                    run_prefix, read_ene=True, read_zpe=True, saddle=False)
            else:
                ene_spc, ene_basis = thmroutines.basis.basis_energy(
                    spc_name, spc_basis, uniref_dct, spc_dct,
                    pf_levels[spc_model], pf_models[spc_model],
                    run_prefix, save_prefix)
                for spc_basis_i, ene_basis_i in zip(spc_basis, ene_basis):
                    chn_basis_ene_dct[spc_model][spc_basis_i] = ene_basis_i

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
            ckin_nasa_str += writer.ckin.model_header(
                spc_models,
                pf_levels[spc_model],
                pf_models[spc_model],
                refscheme=spc_model_dct[spc_model]['therm']['ref_scheme'])

            # Build and write the NASA polynomial in CHEMKIN-format string
            ckin_nasa_str += thmroutines.nasapoly.build_polynomial(
                spc_name, spc_dct,
                thm_paths[idx]['final'][0], thm_paths[idx]['final'][1])
            ckin_nasa_str += '\n\n'

            print(ckin_nasa_str)

        # Write all of the NASA polynomial strings
        writer.ckin.write_nasa_file(ckin_nasa_str, ckin_path)
