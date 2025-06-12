""" test autofile.fs
"""
import os
import tempfile
# import numpy
# import pytest

import automol
import autofile.fs

PREFIX = tempfile.mkdtemp()
print(PREFIX)


def test__manager():
    """ test
    """

    prefix = os.path.join(PREFIX, 'data1')
    _build_fs(prefix)

    # Set locs for this test
    spc_locs = ['InChI=1S/H2O/h1H2', 0, 1]
    thy_locs = ['hf', 'sto-3g', 'U']
    cnf_locs_1 = ['rQ5VxakIXDkDp', 'cdgZx6pwjFtcX']
    cnf_locs_2 = ['rTW2L_4wU-PVV', 'caI7Oi5445ya2']

    # Set cnf fs using the manager
    cnf_fs = autofile.fs.manager(
        prefix, [['SPECIES', spc_locs], ['THEORY', thy_locs]], 'CONFORMER')

    # Test cnf fs with removing files
    # Removing files/layers withut setting removable
    rem1, rem2 = False, False
    try:
        cnf_fs[-1].remove(cnf_locs_1)
    except ValueError:
        rem1 = True
    try:
        cnf_fs[-1].file.geometry_input.remove(cnf_locs_1)
    except ValueError:
        rem2 = True
    assert rem1 and rem2

    cnf_fs[-1].removable = True
    if cnf_fs[-1].file.geometry_input.exists(cnf_locs_2):
        cnf_fs[-1].file.geometry_input.removable = True
        cnf_fs[-1].file.geometry_input.remove(cnf_locs_2)
    assert not cnf_fs[-1].file.geometry_input.exists(cnf_locs_2)
    cnf_fs[-1].remove(cnf_locs_2)


def test__iterate_managers():
    """ test autofile.fs.iterate_managers
    """

    def _get_set(locs):
        """ Figure out set name from cnf locs
        """
        for key, val in FAKE_LOCS_SETS_DCT.items():
            if locs == val[2]:
                set_name = key
                break
        return set_name

    # Build fs
    prefix = os.path.join(PREFIX, 'data2')
    _build_fs(prefix)

    # Generate managers for grapping parts of fake fileystem
    cnf_managers = tuple(autofile.fs.iterate_managers(
        prefix, ['SPECIES', 'THEORY'], 'CONFORMER'))

    # Loop over all cnf_fs and read file for correct data
    for cnf_fs in cnf_managers:
        for locs in cnf_fs[-1].existing():
            inp_str = cnf_fs[-1].file.geometry_input.read(locs)
            assert inp_str == FAKE_NAMES_DCT[_get_set(list(locs))]


def test__iterate_locators():
    """ test autofile.fs.iterate_locators
    """

    # Build fs
    prefix = os.path.join(PREFIX, 'data3')
    _build_fs(prefix)

    # Generate managers for grapping parts of fake fileystem
    spc_locators = tuple(autofile.fs.iterate_locators(
        prefix, ['SPECIES']))

    # Loop over species locs (generate geo from ich and match filesys)
    for spc_locs, in spc_locators:
        ich, _, _ = spc_locs
        ref_geo = automol.chi.geometry(ich)

        spc_fs = autofile.fs.species(prefix)
        thy_fs = autofile.fs.theory(spc_fs[-1].path(spc_locs))
        for thy_locs in thy_fs[-1].existing():
            cnf_fs = autofile.fs.conformer(thy_fs[-1].path(thy_locs))
            for cnf_locs in cnf_fs[-1].existing():
                geo = cnf_fs[-1].file.geometry.read(cnf_locs)
                assert automol.geom.almost_equal_dist_matrix(ref_geo, geo)


def _build_fs(prefix):
    """ Construct a filesystem to do stuff
    """

    os.mkdir(prefix)

    spc_fs = autofile.fs.species(prefix)
    for _set in ('set1', 'set2', 'set3', 'set4'):
        locs_set = FAKE_LOCS_SETS_DCT[_set]
        inp_name = FAKE_NAMES_DCT[_set]
        # spc layer
        spc_fs[-1].create(locs_set[0])
        spc_path = spc_fs[-1].path(locs_set[0])
        # thy layer
        thy_fs = autofile.fs.theory(spc_path)
        thy_fs[-1].create(locs_set[1])
        thy_path = thy_fs[-1].path(locs_set[1])
        # cnf layer
        cnf_fs = autofile.fs.conformer(thy_path)
        cnf_fs[-1].create(locs_set[2])
        cnf_fs[-1].file.geometry_input.write(inp_name, locs_set[2])
        geo = automol.chi.geometry(locs_set[0][0])
        cnf_fs[-1].file.geometry.write(geo, locs_set[2])


FAKE_LOCS_SETS_DCT = {
    'set1': [['InChI=1S/H2O/h1H2', 0, 1],
             ['hf', 'sto-3g', 'U'],
             ['rQ5VxakIXDkDp', 'cdgZx6pwjFtcX']],
    'set2': [['InChI=1S/H2O/h1H2', 0, 1],
             ['hf', 'sto-3g', 'U'],
             ['rTW2L_4wU-PVV', 'caI7Oi5445ya2']],
    'set3': [['InChI=1S/CH4/h1H4', 0, 1],
             ['hf', 'cc-pvdz', 'R'],
             ['r5-8BlffIKsw_', 'cWgRuV5qGsBrK']],
    'set4': [['InChI=1S/CH4/h1H4', 0, 1],
             ['hf', 'cc-pvdz', 'R'],
             ['ry0F_raMFxFnv', 'ckyd_uZJ-TmM7']]
}
FAKE_NAMES_DCT = {
    'set1': 'inp1',
    'set2': 'inp2',
    'set3': 'inp3',
    'set4': 'inp4'
}
