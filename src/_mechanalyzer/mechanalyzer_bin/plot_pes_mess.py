""" Script to read in a PES from a MESS input file and then plot it.

    PLot generated is ReactionCoordinate vs. Relative Energy.
"""

import sys
import os
import argparse
import ioformat
import automol
import mess_io.reader
import mechanalyzer.plotter


# Set useful global variables
CWD = os.getcwd()

# Parse the command line
DSTR = 'Generates a graph plot of a PES from MESS input'
PAR = argparse.ArgumentParser(description=DSTR)
PAR.add_argument('-i', '--input', default='mess.inp',
                 help='input file type (mess.inp')
PAR.add_argument('-o', '--output', default='surface.pdf',
                 help='name of output plot file (surface.pdf)')
OPTS = vars(PAR.parse_args())

# Parse the input based on the initial type
INP_PES_STR = ioformat.pathtools.read_file(CWD, OPTS['input'])
if INP_PES_STR is None:
    print(f'ERROR: Input PES file {OPTS["input"]} not found')
    sys.exit()

ene_dct, _, conn_lst_dct, pes_lab_dct = mess_io.reader.pes(
    INP_PES_STR, read_fake=False)
ene_dct = {name: ene for name, ene in ene_dct.items()
           if 'B' not in name}
conn_lst = tuple(conn[1] for conn in conn_lst_dct.items())
pes_lab_dct = automol.util.dict_.invert(pes_lab_dct)
print(pes_lab_dct)

# Call the plotter function
print('Information parsed from MESS input file')
print('Energies of species:')
for name, ene in ene_dct.items():
    print(f'{name}: {ene} kcal/mol')

print('Connections:')
for conn in conn_lst:
    print(f'{conn}')

# Produce the PES plot
mechanalyzer.plotter.pes.pes_graph(
    conn_lst, ene_dct=ene_dct, label_dct=pes_lab_dct)

# Exit
print('\nPlot surface created. Script Execution complete.')
