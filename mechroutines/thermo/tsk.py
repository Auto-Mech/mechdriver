""" Tasks for THERMODRIVER
"""

import autorun
from automol.chi import formula_string as fstring
import thermfit
from mechlib import filesys
from mechlib.amech_io import reader
from mechlib.amech_io import writer
from mechlib.amech_io import parser
from mechlib.amech_io import output_path
from mechlib.amech_io import printer as ioprinter
from mechroutines.models import ene
from mechroutines.thermo import qt
from mechroutines.thermo import nasapoly
from mechroutines.thermo import basis as thmbasis


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
                messpf_inp_str, dat_dct = qt.make_messpf_str(
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

    ioprinter.info_message('\n\n')
    ioprinter.obj('line_dash')


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
        ioprinter.message(f'Run MESSPF: {spc_name}', newline=1)
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

    ioprinter.info_message('\n\n')
    ioprinter.obj('line_dash')


def produce_boltzmann_weighted_conformers_pf(
        run_messpf_tsk, spc_locs_dct, spc_dct,
        thm_paths_dct):
    """ Combine PFs into final pf
    """
    ioprinter.messpf('run_header')

    spc_mods, _ = parser.models.extract_models(run_messpf_tsk)
    print('starting produce_boltz...')

    for spc_name in spc_locs_dct:
        ioprinter.message(f'Run MESSPF: {spc_name}', newline=1)
        locs_pfs_arrays = []
        hf_array = []
        for idx, spc_locs in enumerate(spc_locs_dct[spc_name]):
            locs_pfs_arrays.append(reader.mess.messpf(
                thm_paths_dct[spc_name][tuple(spc_locs)]['mod_total'][0]))
            hf_val = 0.
            for spc_mod in spc_mods:
                hf_val += (
                    spc_dct[spc_name]['Hfs'][idx][spc_mod][0] / len(spc_mods)
                )
            hf_array.append(hf_val)
        final_pf = thermfit.pf.boltzmann_pf_combination(
            locs_pfs_arrays, hf_array)
        writer.mess.output(
            fstring(spc_dct[spc_name]['inchi']),
            final_pf,
            thm_paths_dct[spc_name]['spc_total'][0],
            filename='pf.dat')
        spc_dct[spc_name]['Hfs']['final'] = [min(hf_array)]
    return spc_dct


def multi_species_pf(
        run_messpf_tsk, spc_locs_dct, spc_dct,
        thm_paths_dct, spc_grp_lst):
    """ Combine PFs into final pf
    """
    ioprinter.messpf('run_header')

    spc_mods, _ = parser.models.extract_models(run_messpf_tsk)
    print('starting produce_boltz...')

    for x in spc_locs_dct:
        print(x)
    print('---')

    for grp_name, grp_lst in spc_grp_lst.items():
        locs_pfs_arrays_lst = []
        hf_array_lst = []
        for spc_name in grp_lst:
            hf_array = []
            locs_pfs_arrays = []
            ioprinter.message(f'Run MESSPF: {spc_name}', newline=1)
            for idx, spc_locs in enumerate(spc_locs_dct[spc_name]):
                locs_pfs_arrays.append(reader.mess.messpf(
                    thm_paths_dct[spc_name][tuple(spc_locs)]['mod_total'][0]))
                hf_val = 0.
                for spc_mod in spc_mods:
                    hf_val += (
                        spc_dct[spc_name]['Hfs'][idx][spc_mod][0] / len(spc_mods)
                    )
                hf_array.append(hf_val)
                hf_array_lst.append(hf_val)
            final_pf = thermfit.pf.boltzmann_pf_combination(
                locs_pfs_arrays, hf_array)
            locs_pfs_arrays_lst.append(final_pf)
            writer.mess.output(
                fstring(spc_dct[spc_name]['inchi']),
                final_pf,
                thm_paths_dct[spc_name]['spc_total'][0],
                filename='pf.dat')
            print('mess output at', thm_paths_dct[spc_name]['spc_total'][0])
            spc_dct[spc_name]['Hfs']['final'] = [min(hf_array)]
        
        # Get the final species named and formula
        init_spc = grp_lst[0]
        spc_dct[grp_name] = spc_dct[init_spc]
        spc_fml = fstring(spc_dct[grp_name]['inchi'])

        # Combine the PFs and H0K for all confs of all species
        # final_combo_pf = thermfit.pf.combine_pfs_additively(
        #     locs_pfs_arrays_lst)
        final_combo_pf = thermfit.pf.stereo_pf_combination(
            locs_pfs_arrays_lst, hf_array_lst)

        # Write MESSPF file for combined PFs
        print('mess output at', thm_paths_dct[grp_name]['spc_group'][0])
        writer.mess.output(
            spc_fml,
            final_combo_pf,
            thm_paths_dct[grp_name]['spc_group'][0],
            filename='pf.dat')

        spc_dct[grp_name]['Hfs']['final'] = [min(hf_array_lst)]
        # Add the H0K to spc_dct for spc
        # Set a new spc_dct entry
        # if len(spc_grp) > 1:
        #     spc_dct[grp_name] = init_spc_dct
        #     spc_dct[grp_name]['Hfs']['final'] = [min(hf_array)]
        #     thm_paths_dct[grp_name]['spc_total'] = (
        #         job_path(
        #             run_prefix, 'MESS', 'PF',
        #             thm_prefix, locs_id=idx),
        #         job_path(
        #             run_prefix, 'THERM', 'NASA',
        #             thm_prefix, locs_id=idx)
        #     )

    return spc_dct


def _weigh_heat_of_formation(hf_array, weights):
    """ weight heat of formation of conformers based on the weights
        determined during the combination of pfs
    """
    hf_val = 0.
    for hf_val_i, weight_i in zip(hf_array, weights):
        hf_val += weight_i * hf_val_i
    return hf_val


def _get_heat_of_formation(
        spc_name, spc_dct, spc_mod, spc_locs,
        spc_mod_dct_i, basis_dct, ref_enes,
        chn_basis_ene_dct, run_prefix, save_prefix):
    """ get the heat of formation at  0 K
    """

    spc_basis, coeff_basis = basis_dct[spc_name]
    # Get the energies for the spc and its basis
    energy_missing, ene_basis = _check_for_reference_energies(
        spc_basis, chn_basis_ene_dct, spc_mod)
    if energy_missing:
        print('uh oh i messed things up')
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct[spc_name], spc_mod_dct_i,
        run_prefix, save_prefix, saddle=False, spc_locs=spc_locs)
    ene_spc = ene.read_energy(
        spc_dct[spc_name], pf_filesystems, spc_mod_dct_i,
        run_prefix, read_ene=True, read_zpe=True, saddle=False)
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
    # Determine info about the basis species used in thermochem calcs
    for spc_mod in spc_mods:
        spc_mod = spc_mods[0]
        spc_mod_dct_i = spc_mod_dct[spc_mod]
        chn_basis_ene_dct[spc_mod] = {}
        basis_dct = thermfit.prepare_basis(
            ref_scheme, spc_dct, (*spc_locs_dct.keys(),))
        uniref_dct = thermfit.unique_basis_species(basis_dct, spc_dct)
        # uniich_lst = [uniref_dct[key]['inchi'] for key in uniref_dct]
        basis_ichs = []
        for parent_spc_name in basis_dct:
            for basis_spc_name in basis_dct[parent_spc_name][0]:
                if basis_spc_name not in basis_ichs:
                    basis_ichs.append(basis_spc_name)
        _, ene_basis = thmbasis.basis_energy(
            None, basis_ichs,
            uniref_dct, spc_dct,
            spc_mod_dct_i,
            run_prefix, save_prefix, read_species=False)
        for spc_basis_i, ene_basis_i in zip(basis_ichs, ene_basis):
            chn_basis_ene_dct[spc_mod][spc_basis_i] = ene_basis_i

    for spc_name in spc_locs_dct:
        for idx, spc_locs in enumerate(spc_locs_dct[spc_name]):
            spc_locs = tuple(spc_locs)
            for spc_mod in spc_mods:
                print("\n+++++++++++++++++++++++++++++++++++++++++++++\n")
                print(" Calculating 0 K Heat-of-Formation for "
                      f"{spc_name} {spc_mod}\n")
                # Take species model and add it to the chn_basis_ene dct
                spc_mod = spc_mods[0]
                spc_mod_dct_i = spc_mod_dct[spc_mod]
                hf0k, chn_basis_ene_dct = _get_heat_of_formation(
                    spc_name, spc_dct, spc_mod, spc_locs,
                    spc_mod_dct_i, basis_dct, ref_enes,
                    chn_basis_ene_dct, run_prefix, save_prefix)
                spc_dct = _add_hf_to_spc_dct(
                    hf0k, spc_dct, spc_name, idx, spc_mod)
    return spc_dct


def nasa_polynomial_task(
        mdriver_path, spc_locs_dct, thm_paths_dct, spc_dct,
        spc_mod_dct, spc_mods, sort_info_lst, ref_scheme,
        spc_grp_dct=None):
    ckin_nasa_str_dct = {}
    ckin_nasa_str_dct[0] = ''
    """ generate the nasa polynomials
    """
    ckin_nasa_str_dct = {}
    ckin_nasa_str_dct[0] = ''
    ckin_path = output_path('CKIN', prefix=mdriver_path)
    for spc_name in spc_locs_dct:
        for idx, spc_locs in enumerate(spc_locs_dct[spc_name], start=1):
            if idx not in ckin_nasa_str_dct:
                ckin_nasa_str_dct[idx] = ''
            spc_locs = tuple(spc_locs)
            ioprinter.nasa('calculate', spc_name)
            ioprinter.message('for: ', spc_locs, ' combined models')
            ckin_nasa_str_dct[idx] += writer.ckin.model_header(
                spc_mods, spc_mod_dct,
                sort_info_lst=sort_info_lst,
                refscheme=ref_scheme)
            ckin_nasa_str_dct[idx] += nasapoly.build_polynomial(
                spc_name, spc_dct,
                thm_paths_dct[spc_name][tuple(spc_locs)]['mod_total'][0],
                thm_paths_dct[spc_name][tuple(spc_locs)]['mod_total'][1],
                spc_locs_idx=idx-1, spc_mod=','.join(spc_mods))
            ckin_nasa_str_dct[idx] += '\n\n'
            # ioprinter.info_message('CKIN NASA STR\n')
            # ioprinter.info_message(ckin_nasa_str_dct[idx])
        ioprinter.message('for combined rid cids:', spc_locs_dct[spc_name])
        ckin_nasa_str_dct[0] += writer.ckin.model_header(
            spc_mods, spc_mod_dct,
            sort_info_lst=sort_info_lst,
            refscheme=ref_scheme)
        ckin_nasa_str_dct[0] += nasapoly.build_polynomial(
            spc_name, spc_dct,
            thm_paths_dct[spc_name]['spc_total'][0],
            thm_paths_dct[spc_name]['spc_total'][1],
            spc_locs_idx='final', spc_mod=','.join(spc_mods))
        ckin_nasa_str_dct[0] += '\n\n'
    for idx in ckin_nasa_str_dct:
        ioprinter.info_message('CKIN NASA STR {}\n'.format(str(idx)))
        ioprinter.info_message(ckin_nasa_str_dct[idx])

    if spc_grp_dct is not None:
        ckin_nasa_str_dct[1000] = ''
        for grp_name in spc_grp_dct:
            ioprinter.message('for combined species:', grp_name)
            ckin_nasa_str_dct[1000] += writer.ckin.model_header(
                spc_mods, spc_mod_dct,
                sort_info_lst=sort_info_lst,
                refscheme=ref_scheme)
            ckin_nasa_str_dct[1000] += nasapoly.build_polynomial(
                grp_name, spc_dct,
                thm_paths_dct[grp_name]['spc_group'][0],
                thm_paths_dct[grp_name]['spc_group'][1],
                spc_locs_idx='final', spc_mod=','.join(spc_mods))
            ckin_nasa_str_dct[1000] += '\n\n'
        ioprinter.info_message('CKIN NASA STR COMBINED SPECIES\n')
        ioprinter.info_message(ckin_nasa_str_dct[1000])

    return ckin_nasa_str_dct, ckin_path
