""" Script for running a comparison of rate constants between mechanisms
"""

import os
import sys
import numpy
import mechanalyzer.calculator.compare as compare
import mechanalyzer.plotter.rates as plot_rates
import mechanalyzer.plotter._util as plot_util
import mechanalyzer.parser.new_spc as spc_parser
import mechanalyzer.parser.ckin_ as ckin_parser
from mechanalyzer.cli import util
from ioformat import pathtools


def main(
        mechs_yaml,
        plot_fname,
        out_txt_fname,
        job_path,
        temps_lst,
        pressures,
        sort_method,
        rev_rates,
        remove_loners):

    labels, mech_files, therm_files, csv_files = util.read_mechs_yaml(
        mechs_yaml, 'rates')
    if temps_lst is None or len(temps_lst) < 2:
        temps_lst = [numpy.linspace(500, 1500, 16)]
    else:
        temps_lst = [numpy.linspace(temps_lst[0], temps_lst[1], 16)]
    if pressures is None or len(pressures) < 1:
        pressures = [1, 10, 100]

    rxn_ktp_dcts = ckin_parser.load_rxn_ktp_dcts(
        mech_files, job_path, temps_lst, pressures)
    spc_therm_dcts = ckin_parser.load_spc_therm_dcts(
        therm_files, job_path, temps_lst[0])  # NOTE: taking first entry
    spc_dcts = spc_parser.load_mech_spc_dcts(csv_files, job_path)

    # Get the algn_rxn_ktp_dct
    temps = temps_lst[0]  # function receives a single Numpy array of temps
    algn_rxn_ktp_dct = compare.get_algn_rxn_ktp_dct(
        rxn_ktp_dcts, spc_therm_dcts, spc_dcts, temps, rev_rates=rev_rates,
        remove_loners=remove_loners)

    # Run the plotter
    figs, _ = plot_rates.build_plots(
        algn_rxn_ktp_dct,
        mech_names=labels,
        ratio_sort=bool(sort_method == 'ratios'))
    plot_util.build_pdf(figs, filename=plot_fname, path=job_path)

    # Write the ordered text file
    FSTR = compare.write_ordered_str(algn_rxn_ktp_dct, dct_type='rxn')
    pathtools.write_file(FSTR, job_path, out_txt_fname)

