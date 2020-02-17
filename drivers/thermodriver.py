""" Driver for thermochemistry evaluations including
    heats-of-formation and NASA polynomials describing
    thermodynamic quantities: Enthalpy, Entropy, Gibbs
"""

import os
import autofile.fs
import routines
from lib.load import model as loadmodel
from lib.load import mechanism as loadmech
from lib.filesystem import inf as finf
from lib.runner import therm as thmrunner
from lib.outpt import chemkin as cout
from lib import printmsg


def run(spc_dct,
        model_dct,
        thy_dct,
        rxn_lst,
        run_inp_dct,
        ref_scheme='basic',
        run_pf=True,
        run_thermo=True):
    """ main driver for thermo run
    """

    # Print the header message for the driver
    printmsg.program_header('thermo')

    # Pull stuff from dcts for now
    save_prefix = run_inp_dct['save_prefix']
    ref_set = 'ANL0'
    # ene_idx = 0

    # Build a list of the species to calculate thermochem for loops below
    spc_queue = loadmech.build_spc_queue(rxn_lst)

    # Determine information about the basis species used in thermochem calcs
    basis_dct, unique_basis_dct, msg = routines.pf.thermo.therm.prepare_refs(
        ref_scheme, spc_dct, spc_queue)

    # Add the references to the spc_dct
    spc_dct.update(unique_basis_dct)

    print('basis dct')
    for key, val in basis_dct.items():
        print(key)
        print(val)

    # Print the info message
    print(msg)

    # Write and Run MESSPF inputs to generate the partition functions
    if run_pf:

        # Write and Run the MESSPF file to get the partition function
        for spc in spc_queue:

            # Unpack spc to get name and model
            spc_name, spc_model = spc

            # Gather PF model and theory level info
            pf_levels = loadmodel.set_es_model_info(
                model_dct[spc_model]['es'], thy_dct)
            pf_model = loadmodel.set_pf_model_info(
                model_dct[spc_model]['pf'])
            freeze_all_tor = model_dct[spc_model]['options']['freeze_all_tors']
            ndim_tors = model_dct[spc_model]['pf']['tors']

            # Get PF input header
            temp_step = 100.
            ntemps = 30
            global_pf_str = routines.pf.messf.pfblock.get_pf_header(
                temp_step, ntemps)

            # Set up the species filesystem
            spc_info = finf.get_spc_info(spc_dct[spc_name])
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs[-1].create(spc_info)
            spc_save_path = spc_save_fs[-1].path(spc_info)

            # Read the ZPVE from the filesystem if not doing tau
            tau_mod = bool(model_dct[spc_model]['pf']['tors'] == 'tau')
            if not tau_mod:
                zpe = routines.pf.messf.get_zero_point_energy(
                    spc, spc_dct[spc_name], pf_levels, pf_model,
                    save_prefix=spc_save_path)
                zpe_str = routines.pf.messf.get_zpe_str(
                    spc_dct[spc_name], zpe)
            else:
                zpe_str = ''

            # Generate the partition function
            spc_str, data_str_dct, _ = routines.pf.messf.blocks.species_block(
                spc=spc,
                spc_dct_i=spc_dct[spc_name],
                spc_info=spc_info,
                spc_model=pf_model,
                pf_levels=pf_levels,
                save_prefix=spc_save_path,
                tors_mod=(ndim_tors, freeze_all_tor)
                )

            # Write the MESSPF input file
            harm_thy_info = pf_levels[2]
            pf_input = routines.pf.messf.pfblock.get_pf_input(
                spc_name, spc_str, global_pf_str, zpe_str)
            pf_path, nasa_path = thmrunner.get_thermo_paths(
                spc_save_path, spc_info, harm_thy_info)
            routines.pf.messf.pfblock.write_pf_input(
                pf_input, data_str_dct, pf_path)

            # Run the MESSPF File
            thmrunner.run_pf(pf_path)

    # Use MESS partition functions to compute thermo quantities
    if run_thermo:

        # Setup the CHEMKIN level comment string
        # ene_str = cout.get_ckin_ene_lvl_str(pf_levels, ene_coeff)

        # Read the high-level energy
        for spc in spc_queue:

            # Unpack spc to get name and model
            spc_name, spc_model = spc

            # Get the basis info for the spc of interest
            spc_basis, coeff_basis = basis_dct[spc]

            # Get the energies for the spc and its basis
            ene_spc = routines.pf.messf.ene.get_fs_ene_zpe(
                spc_dct, spc_name,
                thy_dct, model_dct, spc_model,
                save_prefix, saddle=False,
                ene_coeff=[1.0],
                read_ene=True, read_zpe=True)
            ene_basis = routines.pf.thermo.therm.basis_energy(
                spc_basis, spc_dct,
                thy_dct, model_dct, spc_model, save_prefix,
                ene_coeff=[1.0])

            # Calculate and store the 0 K Enthalpy
            hf0k = routines.pf.thermo.heatform.calc_hform_0k(
                ene_spc, ene_basis, spc_basis, coeff_basis, ref_set=ref_set)
            spc_dct[spc_name]['Hfs'] = [hf0k]

        # Write the NASA polynomials in CHEMKIN format
        chemkin_header_str = cout.run_ckin_header(
            pf_levels, pf_model)
        chemkin_set_str = chemkin_header_str
        for spc in spc_queue:

            # Unpack spc to get name and model
            spc_name, spc_model = spc

            # Set up the paths for running jobs
            harm_thy_info = pf_levels[2]
            pf_path, nasa_path = thmrunner.get_thermo_paths(
                spc_save_path, spc_info, harm_thy_info)
            starting_path = thmrunner.go_to_path(nasa_path)

            # Write and run the thermp file to get Hf0k and ...
            thmrunner.write_thermp_inp(spc_dct[spc_name])
            # not getting spc str, so this isnt working, fix this
            # if spc_dct[spc]['ene'] == 0.0 or spc_dct[spc]['spc_str'] == '':
            #     print('Cannot generate thermo for species',
            #           '{} '.format(spc_dct[spc]['ich']),
            #           'because information is still missing:')
            #     continue
            hf298k = thmrunner.run_thermp(pf_path, nasa_path)
            spc_dct[spc_name]['Hfs'].append(hf298k)

            # Run PAC99 to get a NASA polynomial string in its format
            pac99_poly_str = thmrunner.run_pac(spc_dct[spc_name], nasa_path)

            # Convert the polynomial from PAC99 to CHEMKIN
            chemkin_poly_str = cout.run_ckin_poly(
                spc_name, spc_dct[spc_name], pac99_poly_str)

            # Write a string for a single spc file to a set
            chemkin_spc_str = chemkin_header_str + chemkin_poly_str
            chemkin_set_str += chemkin_poly_str

            # Write the CHEMKIN string to a file
            thmrunner.go_to_path(starting_path)
            ckin_path = thmrunner.prepare_path(starting_path, 'ckin')
            if not os.path.exists(ckin_path):
                os.makedirs(ckin_path)
            cout.write_nasa_file(
                spc_dct[spc_name], ckin_path, nasa_path, chemkin_spc_str)

        # Write the NASA polynomial strings to a CHEMKIN-formatted file
        with open(os.path.join(ckin_path, 'automech.ckin'), 'w') as nasa_file:
            nasa_file.write(chemkin_set_str)
