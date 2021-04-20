""" Tests the input for the sorter
    Sorts inp/mechanism.dat according
    to inp/sort.dat criteria
"""

import os
import subprocess
from drivers import sortdriver
from mechlib.amech_io import parser

PATH = os.getcwd()
PATHINP = os.path.join(PATH, 'inp')
isolate_species, sort_list = parser.mechanism.read_sort_section(PATH)
sortdriver.run(PATHINP, 'species.csv', 'mechanism.dat',
               'mech_sort.dat', sort_list, isolate_species=isolate_species)
