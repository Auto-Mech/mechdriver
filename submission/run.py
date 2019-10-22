import read_dat

PARAMS = read_dat.params('params.dat')

for attr in dir(PARAMS):
    if '__' not in attr:
        print(attr, getattr(PARAMS, attr), type(getattr(PARAMS, attr)))
