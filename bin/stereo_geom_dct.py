import sys

from autofile import io_ as io
import automol
import pprint

geo_file = sys.argv[1]
geo_str = io.read_file(geo_file)
geo = automol.geom.from_string(geo_str)
ich = automol.geom.inchi(geo, stereo=True)
ich_key = automol.inchi.inchi_key(ich)
formula = automol.geom.formula(geo)
gra = automol.geom.graph(geo)
smiles = automol.inchi.smiles(ich)

bad_ich_dct = {
    ich: {
        'inchi': ich,
        'geom': geo,
        'formula': formula,
        'smiles': smiles,
        'graph': gra}
    }

print('Add this dictionary to autochem/automol/inchi/base/_core.py')
#pprint.pprint(bad_ich_dct)
print(bad_ich_dct)
print('\ncopy your conformer from that UHFFF guy into:  ', ich_key)
print('\nchange the dir.yaml to have the right inchi', ich)
