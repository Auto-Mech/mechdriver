""" Useful functions for plotting
"""

import os
import matplotlib.backends.backend_pdf as plt_pdf


def build_pdf(figs, filename='output.pdf', path=None):
    """ Produce a PDF with one reaction per page

        :param figs: list of MatPlotLib figure objects
        :type figs: list [fig1, fig2, ...]
        :param filename: filename for the output pdf
        :type filename: str
        :param path: path for the output pdf; default is the current directory
        :type path: str
    """
    print('Producing PDF...')
    if path is not None:
        filename = os.path.join(path, filename)
    pdf = plt_pdf.PdfPages(filename)
    for fig in figs:
        pdf.savefig(fig)
    pdf.close()
