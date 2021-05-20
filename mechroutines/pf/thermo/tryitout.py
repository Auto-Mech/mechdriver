""" Script for calculating CBH values
"""

import mechroutines.pf.thermo as thermo
from autofile import io_ as io
from automol import geom

# Set reaction info
frm_key = [5, 6]
brk_key = [1, 5]
RXNCLASS = 'h_abstraction'

geostr = io.read_file('geom.xyz')
geo = geom.from_string(geostr)
frm_key = [5, 6]
brk_key = [1, 5]
zma = geom.zmatrix(geo, [frm_key, brk_key])
RXNCLASS = 'h_abstraction'

thermo.heatform.get_cbhzed_ts(zma, rxnclass, frm_key, brk_key)
