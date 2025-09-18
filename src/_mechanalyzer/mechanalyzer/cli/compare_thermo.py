""" Script for running a comparison of thermo properties between mechanisms
"""

import os
import sys
import numpy
import mechanalyzer.calculator.compare as compare
import mechanalyzer.plotter.thermo as plot_thermo
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
        sort_method,
        sort_temp,
        remove_loners,
        print_missing
        ):

    labels, _, therm_files, csv_files = util.read_mechs_yaml(
        mechs_yaml, 'thermo')
    if temps_lst is None or len(temps_lst) < 2:
        temps_lst = [numpy.linspace(500, 1500, 16)]
    else:
        temps_lst = [numpy.linspace(temps_lst[0], temps_lst[1], 16)]

    spc_therm_dcts = ckin_parser.load_spc_therm_dcts(
        therm_files, job_path, temps_lst[0])  # NOTE: taking first entry
    spc_dcts = spc_parser.load_mech_spc_dcts(csv_files, job_path)

    # Get the algn_spc_therm_dct
    temps = temps_lst[0]  # function receives a single Numpy array of temps
    algn_spc_therm_dct = compare.get_algn_spc_therm_dct(
        spc_therm_dcts, spc_dcts, remove_loners=remove_loners)

    # Get the combined mech_spc_dct (for plotting the InChI & SMILES)
    comb_spc_dct = compare.get_mult_comb_mech_spc_dct(spc_dcts)

    # Run the plotter
    figs, sort_algn_spc_therm_dct,_ = plot_thermo.build_plots(
        algn_spc_therm_dct,
        spc_dct=comb_spc_dct,
        mech_names=labels,
        sort_method=sort_method,
        sort_temp=sort_temp)
    plot_util.build_pdf(figs, filename=plot_fname, path=job_path)

    # Write the ordered text file
    fstr = compare.write_ordered_str(
        sort_algn_spc_therm_dct, 
        dct_type='therm', 
        comb_mech_spc_dct=comb_spc_dct,
        print_missing=print_missing)
    pathtools.write_file(fstr, job_path, out_txt_fname)
