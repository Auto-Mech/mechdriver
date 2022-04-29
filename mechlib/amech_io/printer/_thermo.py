"""
THERMO routine prints
"""

import chemkin_io
import mechanalyzer
from phydat import phycon

from mechlib.amech_io.printer._print import info_message


def therm_paths_messpf_write_locations(
        spc_name, spc_locs_lst, spc_mods, thm_paths_dct):
    """ prints out a table with the path that the messpf thermo input is
        written for each conformer and model
    """
    info_message('MESSPF location table:')
    info_message(f'{"species name":<16}'
                 f'{"rid":<16}{"cid":<16}{"model":<16}{"path":<16}')
    info_message(f'{"============":<16}'
                 f'{"=======":<16}{"=======":<16}{"=======":<16}'
                 f'{"============":<16}')
    for spc_locs in spc_locs_lst:
        for spc_mod in spc_mods:
            path = thm_paths_dct[spc_name][tuple(spc_locs)][spc_mod][0]
            info_message(
                f'{spc_name:<16}'
                f'{spc_locs[0]:<16}'
                f'{spc_locs[1]:<16}'
                f'{spc_mod:<16}'
                f'{path:<16}'
            )


def therm_paths_messpf_run_locations(
        spc_name, spc_locs_lst, spc_mods, thm_paths_dct):
    """ prints out a table with the path that the messpf thermo is calculated
        for each conformer and model
    """
    info_message('MESSPF location table:')
    info_message(f'{"species name":<16}'
                 f'{"rid":<16}{"cid":<16}{"model":<16}{"path":<16}')
    info_message(f'{"============":<16}'
                 f'{"=======":<16}{"=======":<16}{"=======":<16}'
                 f'{"============":<16}')
    for spc_locs in spc_locs_lst:
        for spc_mod in spc_mods:
            path1 = thm_paths_dct[spc_name][tuple(spc_locs)][spc_mod][0]
            path2 = thm_paths_dct[spc_name][tuple(spc_locs)]['mod_total'][0]
            info_message(
                f'{spc_name:<16}'
                f'{spc_locs[0]:<16}'
                f'{spc_locs[1]:<16}'
                f'{spc_mod:<16}'
                f'{path1:<16}'
            )
            info_message(
                f'{spc_name:<16}'
                f'{spc_locs[0]:<16}'
                f'{spc_locs[1]:<16}'
                f'{"mod combo":<16}'
                f'{path2:<16}'
            )
    info_message(
        f'{spc_name:<16}'
        f'{"":<16}'
        f'{"spc combo":<16}'
        f'{"mod combo":<16}'
        f'{thm_paths_dct[spc_name]["spc_total"][0]:<16}'
    )


def print_thermo(spc_dct, ckin_nasa_str, spc_locs_dct, spc_locs_idx, spc_mod):
    """ Generate and print thermo properties with the nasa polynomial string
    """

    nasa7_params_all = chemkin_io.parser.thermo.create_spc_nasa7_dct(
        ckin_nasa_str)
    templist = (
        298.15, 300, 400, 500, 600, 700, 800,
        900, 1000, 1100, 1200, 1300, 1400, 1500)
    # templist = (
    #     298.15, 300, 400, 500, 600,  800,
    #     1000, 1500, 2000, 2500, 3000)
    info_message(
        'SPECIES            H0f(0 K)  H0f(298 K) in kcal/mol:')
    info_message(
        ' T (K)   H - H(T)    S(T)      Cp(T) ')
    info_message(
        'Kelvin  kcal/mol cal/(mol K) cal/(mol K)')
    whitespace2 = 45
    whitespace2 = whitespace2*' '
    for spc_name in nasa7_params_all:
        nasa7_params = nasa7_params_all[spc_name]
        whitespace = 18-len(spc_name)
        whitespace = whitespace*' '
        if spc_locs_idx in [0, 1000]:
            hf0 = (
                spc_dct[spc_name]['Hfs']['final'][0]
                * phycon.EH2KCAL)
            info_message(
                f'{spc_name}---{"boltzmann_weighted_combo"}')
        else:
            hf0 = (
                spc_dct[spc_name]['Hfs'][spc_locs_idx-1][spc_mod][0]
                * phycon.EH2KCAL)
            idx_str = '_'.join(spc_locs_dct[spc_name][spc_locs_idx-1])
            info_message(f'{spc_name}---{idx_str}')
        hf298 = mechanalyzer.calculator.thermo.enthalpy(
            nasa7_params, 298.15) / 1000.
        info_message(
            f'{whitespace2}{hf0:>9.2f}{hf298:>9.2f}')
        hincref = hf298
        for temp in templist:
            hinct = mechanalyzer.calculator.thermo.enthalpy(
                nasa7_params, temp) / 1000. - hincref
            entt = mechanalyzer.calculator.thermo.entropy(
                nasa7_params, temp)
            cpt = mechanalyzer.calculator.thermo.heat_capacity(
                nasa7_params, temp)
            info_message(
                f'{temp:>8.2f}{hinct:>9.2f}{entt:>9.2f}{cpt:>9.2f}')
