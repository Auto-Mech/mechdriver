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
    info_message('{:<16}{:<16}{:<16}{:<16}{:<16}'.format(
        'species name', 'rid', 'cid', 'model', 'path'))
    info_message('{:<16}{:<16}{:<16}{:<16}{:<16}'.format(
        '============', '=======', '=======', '=======', '=========='))
    for spc_locs in spc_locs_lst:
        for spc_mod in spc_mods:
            info_message('{:<16}{:<16}{:<16}{:<16}{:<16}'.format(
                spc_name, spc_locs[0], spc_locs[1], spc_mod,
                thm_paths_dct[spc_name][tuple(spc_locs)][spc_mod][0]))


def therm_paths_messpf_run_locations(
        spc_name, spc_locs_lst, spc_mods, thm_paths_dct):
    """ prints out a table with the path that the messpf thermo is calculated
        for each conformer and model
    """
    info_message('MESSPF location table:')
    info_message('{:<16}{:<16}{:<16}{:<16}{:<16}'.format(
        'species name', 'rid', 'cid', 'model', 'path'))
    info_message('{:<16}{:<16}{:<16}{:<16}{:<16}'.format(
        '============', '=======', '=======', '======', '=========='))
    for spc_locs in spc_locs_lst:
        for spc_mod in spc_mods:
            info_message('{:<16}{:<16}{:<16}{:<16}{:<16}'.format(
                spc_name, spc_locs[0], spc_locs[1], spc_mod,
                thm_paths_dct[spc_name][tuple(spc_locs)][spc_mod][0]))
            info_message('{:<16}{:<16}{:<16}{:<16}{:<16}'.format(
                spc_name, spc_locs[0], spc_locs[1], 'mod combo',
                thm_paths_dct[spc_name][tuple(spc_locs)]['mod_total'][0]))
    info_message('{:<16}{:<16}{:<16}{:<16}{:<16}'.format(
        spc_name, 'spc combo', '', 'mod combo',
        thm_paths_dct[spc_name]['spc_total'][0]))


def print_thermo(spc_dct, ckin_nasa_str, spc_locs_dct, spc_mod):
    """ Generate and print thermo properties with the nasa polynomial string
    """

    spc_label_dct = {}
    for spc_name in spc_locs_dct:
        for spc_locs in spc_locs_dct[spc_name]:
            spc_label_dct[spc_name + '_' + spc_locs[1][:5]] = [
                spc_name, tuple(spc_locs)]
    nasa7_params_all = chemkin_io.parser.thermo.create_spc_nasa7_dct(
        ckin_nasa_str)
    templist = (
        298.15, 300, 400, 500, 600, 700, 800,
        900, 1000, 1100, 1200, 1300, 1400, 1500)
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
        spc_name_actual, spc_locs = spc_label_dct[spc_name]
        hf0 = (
            spc_dct[spc_name_actual]['Hfs'][spc_locs][spc_mod][0]
            * phycon.EH2KCAL)
        hf298 = mechanalyzer.calculator.thermo.enthalpy(
            nasa7_params, 298.15) / 1000.
        info_message(
            '{}---{}'.format(spc_name_actual, '_'.join(spc_locs)))
        info_message(
            '{}{:>9.2f}{:>9.2f}'.format(whitespace2, hf0, hf298))
        hincref = hf298
        for temp in templist:
            hinct = mechanalyzer.calculator.thermo.enthalpy(
                nasa7_params, temp) / 1000. - hincref
            entt = mechanalyzer.calculator.thermo.entropy(
                nasa7_params, temp)
            cpt = mechanalyzer.calculator.thermo.heat_capacity(
                nasa7_params, temp)
            info_message(
                '{:>8.2f}{:>9.2f}{:>9.2f}{:>9.2f}'
                .format(temp, hinct, entt, cpt))
