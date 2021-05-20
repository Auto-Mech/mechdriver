""" Script for calculating CBH values
"""

import mechroutines.pf.thermo as thermo
from autofile import io_ as io
from automol import geom

# Set Reaction Info
frm_key = [5, 6]
brk_key = [1, 5]
RXNCLASS = 'h_abstraction'

# Get the Z-Matrix
geostr = io.read_file('geom.xyz')
geo = geom.from_string(geostr)
zma = geom.zmatrix(geo, [frm_key, brk_key])

thermo.heatform.get_cbhzed_ts(zma, RXNCLASS, frm_key, brk_key)
