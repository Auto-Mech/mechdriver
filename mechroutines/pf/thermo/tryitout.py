import routines.pf.thermo as thermo
from autofile import io_ as io
from automol import geom


geostr = io.read_file('geom.xyz')
geo = geom.from_string(geostr)
zma = geom.zmatrix(geo, [[1,5],[5,6]])
frm_key = [5, 6]
brk_key = [1, 5]
rxnclass = 'h_abstraction'

thermo.heatform.get_cbhzed_ts(zma, rxnclass, frm_key, brk_key) 

