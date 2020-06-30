""" Driver for thermochemistry evaluations including
    heats-of-formation and NASA polynomials describing
    thermodynamic quantities: Enthalpy, Entropy, Gibbs
"""

import os
from routines.pf import thermo as thmroutines
from routines.pf import runner as pfrunner
from lib.amech_io import parser
from lib.amech_io import writer


def run(spc_dct,
        pes_model_dct, spc_model_dct,
        thy_dct,
        rxn_lst,
        run_inp_dct,
        write_messpf=True,
        run_messpf=True,
        run_nasa=True):
    """ main driver for thermo run
    """

    # Pull stuff from dcts for now
    save_prefix = run_inp_dct['save_prefix']
    run_prefix = run_inp_dct['run_prefix']

    # Build a list of the species to calculate thermochem for loops below
    spc_queue = parser.species.build_queue(rxn_lst)

    # Build the paths [(messpf, nasa)], models and levels for each spc
    starting_path = os.getcwd()
    thm_paths = [pfrunner.thermo_paths(spc_dct[name], run_prefix)
                 for name, _ in spc_queue]
    pf_levels = [parser.pf_level_info(spc_model_dct[spc_model]['es'], thy_dct)
                 for _, (_, spc_model) in spc_queue]
    pf_models = [parser.pf_model_info(spc_model_dct[spc_model]['pf'])
                 for _, (_, spc_model) in spc_queue]

    # Write and Run MESSPF inputs to generate the partition functions
    if write_messpf:

        for idx, (spc_name, (pes_model, spc_model)) in enumerate(spc_queue):

            print("Preparing MESSPF input for ", spc_name)

            global_pf_str = thmroutines.qt.pf_header(
                pes_model_dct[pes_model]['therm_temps'])
            spc_str, dat_str_dct = thmroutines.qt.make_spc_mess_str(
                spc_dct[spc_name], spc_name,
                pf_models[idx], pf_levels[idx],
                run_prefix, save_prefix)
            messpf_inp_str = thmroutines.qt.make_messpf_str(
                global_pf_str, spc_str)
            pfrunner.mess.write_mess_file(
                messpf_inp_str, dat_str_dct, thm_paths[idx][0],
                fname='mess.inp', overwrite=True)

    # Run the MESSPF files that have been written
    if run_messpf:

        for idx, (spc_name, _) in enumerate(spc_queue):

            print("Starting MESSPF calculation for ", spc_name)

            pfrunner.run_pf(thm_paths[idx][0])

    # Use MESS partition functions to compute thermo quantities
    if run_nasa:

        for idx, (spc_name, (pes_model, spc_model)) in enumerate(spc_queue):

            print("Starting thermo calculation for ", spc_name)

            # Get the reference scheme and energies
            ref_scheme = spc_model_dct[spc_model]['options']['ref_scheme']
            ref_enes = spc_model_dct[spc_model]['options']['ref_enes']

            # Determine info about the basis species used in thermochem calcs
            basis_dct, uniref_dct, msg = thmroutines.basis.prepare_refs(
                ref_scheme, spc_dct, spc_queue)
            print(msg)

            # Get the basis info for the spc of interest
            spc_basis, coeff_basis = basis_dct[spc_name]

            # Get the energies for the spc and its basis
            ene_spc, ene_basis = thmroutines.basis.basis_energy(
                spc_name, spc_basis, uniref_dct, spc_dct,
                pf_levels[idx], pf_models[idx],
                run_prefix, save_prefix)

            # Calculate and store the 0 K Enthalpy
            hf0k = thmroutines.heatform.calc_hform_0k(
                ene_spc, ene_basis, spc_basis, coeff_basis, ref_set=ref_enes)
            spc_dct[spc_name]['Hfs'] = [hf0k]

        # Write the NASA polynomials in CHEMKIN format
        ckin_nasa_str = ''
        ckin_path = pfrunner.ckin.path()
        for idx, (spc_name, (pes_model, spc_model)) in enumerate(spc_queue):

            print("Starting NASA polynomials calculation for ", spc_name)

            # Read the temperatures from the pf.dat file, check if viable
            temps = pfrunner.read_messpf_temps(thm_paths[idx][0])
            thmroutines.print_nasa_temps(temps)

            # Write the NASA polynomial in CHEMKIN-format string
            ckin_nasa_str += writer.ckin.model_header(
                pf_levels[idx], pf_models[idx])

            # Build POLY
            ckin_nasa_str += thmroutines.nasapoly.build_polynomial(
                spc_name, spc_dct, temps,
                thm_paths[idx][0], thm_paths[idx][1], starting_path)

        # Write all of the NASA polynomial strings
        writer.chemkin.write_nasa_file(ckin_path, ckin_nasa_str)
