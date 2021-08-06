""" Tasks for THERMODRIVER
"""

import autorun
from automol.inchi import formula_string as fstring
import thermfit
from mechroutines import thermo as thmroutines
from mechroutines.models import ene
from mechlib import filesys
from mechlib.amech_io import reader
from mechlib.amech_io import writer
from mechlib.amech_io import parser
from mechlib.amech_io import output_path
from mechlib.amech_io import printer as ioprinter


def write_messpf_task(
        write_messpf_tsk, spc_locs_dct, spc_dct,
        pes_mod_dct, spc_mod_dct,
        run_prefix, save_prefix, thm_paths_dct):
    """ Write messpf input file
    """
    ioprinter.messpf('write_header')

    spc_mods, pes_mod = parser.models.extract_models(write_messpf_tsk)

    for spc_name in spc_locs_dct:
        ioprinter.therm_paths_messpf_write_locations(
            spc_name, spc_locs_dct[spc_name], spc_mods, thm_paths_dct)
        for spc_locs in spc_locs_dct[spc_name]:
            for spc_mod in spc_mods:
                messpf_inp_str, dat_dct = thmroutines.qt.make_messpf_str(
                    pes_mod_dct[pes_mod]['therm_temps'],
                    spc_dct, spc_name, spc_locs,
                    pes_mod_dct[pes_mod], spc_mod_dct[spc_mod],
                    run_prefix, save_prefix)
                ioprinter.messpf('input_string')
                ioprinter.info_message(messpf_inp_str)
                autorun.write_input(
                    thm_paths_dct[spc_name][tuple(spc_locs)][spc_mod][0],
                    messpf_inp_str,
                    aux_dct=dat_dct,
                    input_name='pf.inp')


def run_messpf_task(
        run_messpf_tsk, spc_locs_dct, spc_dct,
        thm_paths_dct):
    """ Run messpf input file
    """
    ioprinter.messpf('run_header')

    spc_mods, _ = parser.models.extract_models(run_messpf_tsk)

    for spc_name in spc_locs_dct:
        ioprinter.therm_paths_messpf_run_locations(
            spc_name, spc_locs_dct[spc_name], spc_mods, thm_paths_dct)
        # Run MESSPF for all requested models, combine the PFS at the end
        ioprinter.message('Run MESSPF: {}'.format(spc_name), newline=1)
        _locs_pfs = []
        for spc_locs in spc_locs_dct[spc_name]:
            _mod_pfs = []
            for spc_mod in spc_mods:
                autorun.run_script(
                    autorun.SCRIPT_DCT['messpf'],
                    thm_paths_dct[spc_name][tuple(spc_locs)][spc_mod][0])
                _mod_pfs.append(
                    reader.mess.messpf(
                        thm_paths_dct[spc_name][tuple(spc_locs)][spc_mod][0]))

            # Unpack the the pf model combination information
            spc_mod_info = parser.models.split_model(spc_mod)
            _spc_mods, coeffs, operators = spc_mod_info

            final_pf = thermfit.pf.combine(_mod_pfs, coeffs, operators)
            writer.mess.output(
                fstring(spc_dct[spc_name]['inchi']),
                final_pf,
                thm_paths_dct[spc_name][tuple(spc_locs)]['mod_total'][0],
                filename='pf.dat')
            _locs_pfs.append(final_pf)
        # locs_coeffs = [1./len(spc_locs)]*len(spc_locs)
        # locs_operators = ['divide']*(len(spc_locs)-1)
        # final_pf = thermfit.pf.combine(_locs_pfs, locs_coeffs, locs_operators)
        # writer.mess.output(
        #     fstring(spc_dct[spc_name]['inchi']),
        #     final_pf,
        #     thm_paths_dct[spc_name]['spc_total'][0],
        #     filename='pf.dat')


def _get_heat_of_formation(
        spc_name, spc_dct, spc_mod, spc_locs,
        spc_mod_dct_i, ref_scheme, ref_enes,
        chn_basis_ene_dct, run_prefix, save_prefix):
    """ get the heat of formation at  0 K
    """

    # Determine info about the basis species used in thermochem calcs
    basis_dct = thermfit.prepare_basis(
        ref_scheme, spc_dct, (spc_name,))
    uniref_dct = thermfit.unique_basis_species(basis_dct, spc_dct)

    # Get the basis info for the spc of interest
    spc_basis, coeff_basis = basis_dct[spc_name]

    # Get the energies for the spc and its basis
    energy_missing, ene_basis = _check_for_reference_energies(
        spc_basis, chn_basis_ene_dct, spc_mod)

    # if the energies aren't missing we only read the energy
    # of the target species, otherwise we run a function
    # from thermroutines that gathers it for the references
    # as well
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct[spc_name], spc_mod_dct_i,
        run_prefix, save_prefix, saddle=False, spc_locs=spc_locs)
    ene_spc = ene.read_energy(
        spc_dct[spc_name], pf_filesystems, spc_mod_dct_i,
        run_prefix, read_ene=True, read_zpe=True, saddle=False)
    if energy_missing:
        _, ene_basis = thmroutines.basis.basis_energy(
            spc_name, spc_basis, uniref_dct, spc_dct,
            spc_mod_dct_i,
            run_prefix, save_prefix, read_species=False)
        for spc_basis_i, ene_basis_i in zip(spc_basis, ene_basis):
            chn_basis_ene_dct[spc_mod][spc_basis_i] = ene_basis_i

    # Calculate and store the 0 K Enthalpy
    hf0k = thermfit.heatform.calc_hform_0k(
        ene_spc, ene_basis, spc_basis, coeff_basis, ref_set=ref_enes)
    return hf0k, chn_basis_ene_dct


def _check_for_reference_energies(spc_basis, chn_basis_ene_dct, spc_mod):
    """ check the chn_basis_ene_dct for each
        reference species for this spc model
        and return them or say we need to
        prepare references
    """
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
    return energy_missing, ene_basis


def _add_hf_to_spc_dct(hf0k, spc_dct, spc_name, spc_locs_idx, spc_mod):
    """ Put a new hf0k in the species dictioanry
    """
    if 'Hfs' not in spc_dct[spc_name]:
        spc_dct[spc_name]['Hfs'] = {}
    if spc_locs_idx not in spc_dct[spc_name]['Hfs']:
        spc_dct[spc_name]['Hfs'][spc_locs_idx] = {}
    spc_dct[spc_name]['Hfs'][spc_locs_idx][spc_mod] = [hf0k]
    return spc_dct


def get_heats_of_formation(
        spc_locs_dct, spc_dct, spc_mods, spc_mod_dct,
        ref_scheme, ref_enes, run_prefix, save_prefix):
    """ gather the hof
    """
    chn_basis_ene_dct = {}
    for spc_name in spc_locs_dct:
        for idx, spc_locs in enumerate(spc_locs_dct[spc_name]):
            spc_locs = tuple(spc_locs)
            for spc_mod in spc_mods:
                # Take species model and add it to the chn_basis_ene dct
                spc_mod = spc_mods[0]
                spc_mod_dct_i = spc_mod_dct[spc_mod]
                if spc_mod not in chn_basis_ene_dct:
                    chn_basis_ene_dct[spc_mod] = {}
                hf0k, chn_basis_ene_dct = _get_heat_of_formation(
                    spc_name, spc_dct, spc_mod, spc_locs,
                    spc_mod_dct_i, ref_scheme, ref_enes,
                    chn_basis_ene_dct, run_prefix, save_prefix)
                spc_dct = _add_hf_to_spc_dct(
                    hf0k, spc_dct, spc_name, idx, spc_mod)
                # TODO Make a combined Hf0K boltmann or something
    return spc_dct


def nasa_polynomial_task(
        mdriver_path, spc_locs_dct, thm_paths_dct, spc_dct,
        spc_mod_dct, spc_mods, ref_scheme):
    """ generate the nasa polynomials
    """
    ckin_nasa_str_dct = {}
    ckin_path = output_path('CKIN', prefix=mdriver_path)
    for spc_name in spc_locs_dct:
        for idx, spc_locs in enumerate(spc_locs_dct[spc_name]):
            if idx not in ckin_nasa_str_dct:
                ckin_nasa_str_dct[idx] = ''
            spc_locs = tuple(spc_locs)
            ioprinter.nasa('calculate', spc_name)
            # for spc_mod in spc_mods:
            #     ioprinter.message('for: ', spc_locs, spc_mod)
            #     # Write the header describing the models used in thermo calcs
            #     ckin_nasa_str += writer.ckin.model_header(
            #         [spc_mod], spc_mod_dct, refscheme=ref_scheme)
            #     # Build and write the NASA polynomial in CHEMKIN-format string
            #     # Call dies if you haven't run "write mess" task
            #     ckin_nasa_str += thmroutines.nasapoly.build_polynomial(
            #         spc_name, spc_dct,
            #         thm_paths_dct[spc_name][tuple(spc_locs)][spc_mod][0],
            #         thm_paths_dct[spc_name][tuple(spc_locs)][spc_mod][1],
            #         spc_locs=spc_locs, spc_mod=spc_mod)
            #     ckin_nasa_str += '\n\n'
            ioprinter.message('for: ', spc_locs, ' combined models')
            ckin_nasa_str_dct[idx] += writer.ckin.model_header(
                spc_mods, spc_mod_dct, refscheme=ref_scheme)
            ckin_nasa_str_dct[idx] += thmroutines.nasapoly.build_polynomial(
                spc_name, spc_dct,
                thm_paths_dct[spc_name][tuple(spc_locs)]['mod_total'][0],
                thm_paths_dct[spc_name][tuple(spc_locs)]['mod_total'][1],
                spc_locs_idx=idx, spc_mod=','.join(spc_mods))
            ckin_nasa_str_dct[idx] += '\n\n'
        # ioprinter.message('for combined rid cids:', spc_locs_dct[spc_name])
        # ckin_nasa_str += writer.ckin.model_header(
        #     spc_mods, spc_mod_dct, refscheme=ref_scheme)
        # ckin_nasa_str += thmroutines.nasapoly.build_polynomial(
        #     spc_name, spc_dct,
        #     thm_paths_dct[spc_name]['spc_total'][0],
        #     thm_paths_dct[spc_name]['spc_total'][1],
        #     spc_locs=tuple(spc_locs_dct[spc_name][0]), spc_mod=spc_mods[0])
        # ckin_nasa_str += '\n\n'
            ioprinter.info_message('CKIN NASA STR\n')
            ioprinter.info_message(ckin_nasa_str_dct[idx])
    return ckin_nasa_str_dct, ckin_path
