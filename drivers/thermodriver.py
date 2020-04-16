""" Driver for thermochemistry evaluations including
    heats-of-formation and NASA polynomials describing
    thermodynamic quantities: Enthalpy, Entropy, Gibbs
"""

import os
import numpy
import routines
import autofile
from routines.pf import thermo
from lib.load import model as loadmodel
from lib.load import mechanism as loadmech
from lib.filesystem import inf as finf
from lib.runner import therm as thmrunner
from lib.outpt import chemkin as cout
from lib import printmsg


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

    # Print the header message for the driver
    printmsg.program_header('thermo')

    # Pull stuff from dcts for now
    save_prefix = run_inp_dct['save_prefix']

    # Build a list of the species to calculate thermochem for loops below
    spc_queue = loadmech.build_spc_queue(rxn_lst)

    # Write and Run MESSPF inputs to generate the partition functions
    if write_messpf:

        # Write and Run the MESSPF file to get the partition function
        for spc in spc_queue:

            # Unpack spc to get name and model
            spc_name, models = spc
            pes_model, spc_model = models
            print("Preparing messpf input for ", spc_name)

            # Gather PF model and theory level info
            pf_levels = loadmodel.set_es_model_info(
                spc_model_dct[spc_model]['es'], thy_dct)
            pf_model = loadmodel.set_pf_model_info(
                spc_model_dct[spc_model]['pf'])

            # Get PF input header
            temps = pes_model_dct[pes_model]['temps']
            global_pf_str = routines.pf.messf.blocks.get_pf_header(temps)

            # Set up the species filesystem
            spc_info = finf.get_spc_info(spc_dct[spc_name])
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs[-1].create(spc_info)
            spc_save_path = spc_save_fs[-1].path(spc_info)

            # Read the ZPVE from the filesystem if not doing tau
            tau_mod = bool(spc_model_dct[spc_model]['pf']['tors'] == 'tau')
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
                )

            # Write the MESSPF input file
            harm_thy_info = pf_levels[2]
            pf_input = routines.pf.messf.blocks.get_pf_input(
                spc_name, spc_str, global_pf_str, zpe_str)
            pf_path, nasa_path = thmrunner.get_thermo_paths(
                spc_save_path, spc_info, harm_thy_info)
            routines.pf.messf.blocks.write_pf_input(
                pf_input, data_str_dct, pf_path)

    # Use MESS partition functions to compute thermo quantities
    if run_messpf:

        # Setup the CHEMKIN level comment string
        # ene_str = cout.get_ckin_ene_lvl_str(pf_levels, ene_coeff)

        # Read the high-level energy
        for spc in spc_queue:

            # Unpack spc to get name and model
            spc_name, models = spc
            pes_model, spc_model = models
            print("Starting messpf calculation for ", spc_name)

            # Set up the species filesystem
            spc_info = finf.get_spc_info(spc_dct[spc_name])
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs[-1].create(spc_info)
            spc_save_path = spc_save_fs[-1].path(spc_info)

            # Read the ZPVE from the filesystem if not doing tau
            tau_mod = bool(spc_model_dct[spc_model]['pf']['tors'] == 'tau')
            if not tau_mod:
                zpe = routines.pf.messf.get_zero_point_energy(
                    spc, spc_dct[spc_name], pf_levels, pf_model,
                    save_prefix=spc_save_path)
                zpe_str = routines.pf.messf.get_zpe_str(
                    spc_dct[spc_name], zpe)
            else:
                zpe_str = ''

            # Write the MESSPF input file
            harm_thy_info = pf_levels[2]
            pf_path, nasa_path = thmrunner.get_thermo_paths(
                spc_save_path, spc_info, harm_thy_info)

            # Run MESSPF
            thmrunner.run_pf(pf_path)

    # Use MESS partition functions to compute thermo quantities
    if run_nasa:

        # Setup the CHEMKIN level comment string
        # ene_str = cout.get_ckin_ene_lvl_str(pf_levels, ene_coeff)

        # Read the high-level energy
        for spc in spc_queue:

            # Unpack spc to get name and model
            spc_name, models = spc
            pes_model, spc_model = models
            print("Starting thermo calculation for ", spc_name)

            # Get the reference scheme and energies
            ref_scheme = spc_model_dct[spc_model]['options']['ref_scheme']
            ref_enes = spc_model_dct[spc_model]['options']['ref_enes']

            # Determine info about the basis species used in thermochem calcs
            basis_dct, uniref_dct, msg = routines.pf.thermo.basis.prepare_refs(
                ref_scheme, spc_dct, spc_queue)
            print(msg)

            # Get the basis info for the spc of interest
            spc_basis, coeff_basis = basis_dct[spc]

            # Get the energies for the spc and its basis
            ene_spc = routines.pf.messf.ene.get_fs_ene_zpe(
                spc_dct, spc_name,
                thy_dct, spc_model_dct, spc_model,
                save_prefix, saddle=False,
                read_ene=True, read_zpe=True)
            ene_basis = routines.pf.thermo.basis.basis_energy(
                spc_basis, uniref_dct, spc_dct,
                thy_dct, spc_model_dct, spc_model, save_prefix)

            # Calculate and store the 0 K Enthalpy
            hf0k = routines.pf.thermo.heatform.calc_hform_0k(
                ene_spc, ene_basis, spc_basis, coeff_basis, ref_set=ref_enes)
            spc_dct[spc_name]['Hfs'] = [hf0k]

        # Write the NASA polynomials in CHEMKIN format
        chemkin_set_str = ''
        starting_path = thmrunner.get_starting_path()
        ckin_path = thmrunner.prepare_path(starting_path, 'ckin')
        if not os.path.exists(ckin_path):
            os.makedirs(ckin_path)
        for spc in spc_queue:

            # Unpack spc to get name and model
            spc_name, models = spc
            pes_model, spc_model = models
            print("Starting NASA polynomials calculation for ", spc_name)

            # Gather PF model and theory level info
            pf_levels = loadmodel.set_es_model_info(
                spc_model_dct[spc_model]['es'], thy_dct)
            pf_model = loadmodel.set_pf_model_info(
                spc_model_dct[spc_model]['pf'])

            # Begin chemkin string
            # chemkin_header_str = cout.run_ckin_header(
            #     pf_levels, pf_model)
            chemkin_header_str = ''
            chemkin_set_str += chemkin_header_str

            # Set up the paths for running jobs
            harm_thy_info = pf_levels[2]
            pf_path, nasa_path = thmrunner.get_thermo_paths(
                spc_save_path, spc_info, harm_thy_info)
            
            # Read the temperatures from the pf.dat file, check if viable
            temps = thmrunner.read_messpf_temps(pf_path)
            print('Attempting to fit NASA polynomials from',
                    '200-1000 and 1000-3000 K ranges using\n',
                    'Temps from MESSPF file = {}.'.format(
                        (temp for temp in temps
                         if not numpy.isclose(temp, 298.2)))

            # Write and run the thermp file to get Hf0k and ...
            thmrunner.go_to_path(nasa_path)
            thmrunner.write_thermp_inp(spc_dct[spc_name], temps)
            # not getting spc str, so this isnt working, fix this
            # if spc_dct[spc]['ene'] == 0.0 or spc_dct[spc]['spc_str'] == '':
            #     print('Cannot generate thermo for species',
            #           '{} '.format(spc_dct[spc]['ich']),
            #           'because information is still missing:')
            #     continue
            hf298k = thmrunner.run_thermp(pf_path, nasa_path)
            spc_dct[spc_name]['Hfs'].append(hf298k)

            # Run PAC99 to get a NASA polynomial string in its format
            pac99_str = thmrunner.run_pac(spc_dct[spc_name], nasa_path)
            print('str\n', pac99_str)
            pac99_poly_str = thermo.nasapoly.get_pac99_polynomial(pac99_str)

            # Convert the polynomial from PAC99 to CHEMKIN
            chemkin_poly_str = cout.run_ckin_poly(
                spc_name, spc_dct[spc_name], pac99_poly_str)

            # Write a string for a single spc file to a set
            chemkin_spc_str = chemkin_header_str + chemkin_poly_str
            chemkin_set_str += chemkin_poly_str

            # Write the CHEMKIN string to a file
            thmrunner.go_to_path(starting_path)
            cout.write_nasa_file(
                spc_dct[spc_name], ckin_path, nasa_path, chemkin_spc_str)

        # Write the NASA polynomial strings to a CHEMKIN-formatted file
        with open(os.path.join(ckin_path, 'automech.ckin'), 'w') as nasa_file:
            nasa_file.write(chemkin_set_str)
