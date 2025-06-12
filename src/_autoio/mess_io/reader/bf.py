"""
    Read temperature and pressure dependent branching fractions from file
    format as dct
"""

import sys
import numpy as np
import pandas as pd
import copy
import autoparse.find as apf
import autoparse.pattern as app
from ioformat import remove_comment_lines

def get_bf(bf_str):
    
    lines = bf_str.splitlines()
    if '' in lines:
        lines.remove('') #delete empty lines if present
        
    pressures = np.array(lines[0].split()[1:], dtype=float)
    temperatures = np.array([line.split()[0] for line in lines[1:]], dtype=float)
    
    bf_dct = dict.fromkeys(pressures)
    
    for i, p in enumerate(pressures):
        bf = np.array([line.split()[i+1] for line in lines[1:]], dtype=float)
        bf_dct[p] = (temperatures, bf)
        
    return bf_dct
    
