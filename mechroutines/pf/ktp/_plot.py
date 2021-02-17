"""
  Generate a PDF of the plot of a PES
"""

import mess_io
import mechanalyzer
from mechlib.amech_io import printer as ioprinter


def plot_from_dct(ene_dct, conn_lst, formula):
    """ Generate a PES plot
    """
    ioprinter.debug_message('ene', ene_dct)
    ioprinter.debug_message('conn', conn_lst)
    mechanalyzer.plotter.pes.build(ene_dct, conn_lst, formula)


def plot_from_str(mess_input_string, formula):
    """ Generate a PES plot
    """

    # Read the MESS input string and build info needed for plotting
    ene_dct, conn_lst, _ = mess_io.reader.pes(
        mess_input_string, read_fake=False)

    # Generate the plot
    mechanalyzer.plotter.pes.build(ene_dct, conn_lst, formula)
