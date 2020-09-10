"""
  Generate a PDF of the plot of a PES
"""

import mess_io
import mechanalyzer


def plot_pes(mess_input_string, formula):
    """ Generate a PES plot
    """

    # Read the MESS input string and build info needed for plotting
    ene_dct, conn_lst = mess_io.reader.pes(mess_input_string, read_fake=False)

    # Generate the plot
    mechanalyzer.plotter.pes.build(ene_dct, conn_lst, formula)
