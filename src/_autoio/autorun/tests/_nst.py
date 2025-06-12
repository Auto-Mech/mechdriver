""" test autorun.nst
"""

import os
import tempfile
import numpy
import automol
from ioformat import pathtools
import autorun


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')


# Common structural info
CHARGE = 0
MULTS = (1, 3)

# NNO - MOLPRO
NNO_GEO = automol.geom.from_string(
    pathtools.read_file(DAT_PATH, 'nno_init.xyz'))
NNO_PROG = 'molpro2015'
NNO_METHOD = 'caspt2'
NNO_BASIS = 'aug-cc-pvdz'
NNO_ORB_TYPE = 'RR'
NNO_INI_KWARGS = (
    {'gen_lines': {1: ('! wavefunction template',)},
     'casscf_options': ('maxit=40', 'closed,5', 'wf,22,1,0,0'),
     'corr_options': ('shift=0.20',)},
    # 'mol_options': ('nosym',)
    {'gen_lines': {1: ('! wavefunction template',)},
     'casscf_options': ('maxit=40', 'closed,5', 'wf,22,1,2,0'),
     'corr_options': ('shift=0.20',)},
    # 'mol_options': ('nosym',)
)
# ^ need two sets of arguments? have two spin states


# NCN - G09
NCN_GEO = automol.geom.from_string(
    pathtools.read_file(DAT_PATH, 'ncn_init.xyz'))

NCN_PROG = 'gaussian09'
NCN_METHOD = 'mp2'
NCN_BASIS = 'cc-pvdz'
NCN_ORB_TYPE = 'RU'
NCN_INI_KWARGS = (
    {'gen_lines': {1: ('# int=ultrafine',)},
     'machine_options': ('%NProcShared=10',)},
    {'gen_lines': {1: ('# int=ultrafine',)},
     'machine_options': ('%NProcShared=10',)},
)
ZERO = -147.0596823175656


def test__molpro_flux():
    """ test autorun.nst
    """

    run_dir = 'tmp-mol'
    # with tempfile.TemporaryDirectory(dir=PATH) as run_dir:
    isc = autorun.nst.isc_flux(
        run_dir,
        NNO_PROG, NNO_GEO, CHARGE, MULTS,
        NNO_METHOD, NNO_BASIS, NNO_ORB_TYPE, NNO_INI_KWARGS)
    print(isc)


def test__gaussian_flux():
    """ test autorun.nst
    """

    ref_ncn_geo = automol.geom.from_string(
        pathtools.read_file(DAT_PATH, 'ncn_opt.xyz'))
    ref_ncn_flux = pathtools.read_file(DAT_PATH, 'ncn_lz_flux.dat')

    with tempfile.TemporaryDirectory(dir=PATH) as run_dir:
        msx_geo, _, flux_str = autorun.nst.isc_flux(
            run_dir,
            NCN_PROG, NCN_GEO, CHARGE, MULTS,
            NCN_METHOD, NCN_BASIS, NCN_ORB_TYPE, NCN_INI_KWARGS)

        assert automol.geom.almost_equal_dist_matrix(ref_ncn_geo, msx_geo)
        _check_flux_strings(ref_ncn_flux, flux_str)


def _check_flux_strings(ref_flux_str, flux_str):
    """ check the flux strings
    """

    ref_lines = [line.strip() for line in ref_flux_str.splitlines()
                 if line.strip() != '']
    lines = [line.strip() for line in flux_str.splitlines()
             if line.strip() != '']

    for rline, line in zip(ref_lines, lines):
        rene, rflux = rline.split()
        ene, flux = line.split()

        assert numpy.isclose(float(rene), float(ene))
        assert numpy.isclose(float(rflux), float(flux))


if __name__ == '__main__':
    test__molpro_flux()
    # test__gaussian_flux()
