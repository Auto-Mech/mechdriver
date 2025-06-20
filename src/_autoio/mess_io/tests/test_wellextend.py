""" test mess_io._wellextend
"""

import os
import numpy
import automol.util.dict_
from ioformat import pathtools
import mess_io.reader


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH = os.path.join(PATH,'data', 'inp')
OUT_PATH = os.path.join(PATH, 'data', 'out')

KTP_INP_STR = pathtools.read_file(INP_PATH, 'dmm_nowext.inp')
KTP_OUT_STR = pathtools.read_file(OUT_PATH, 'dmm_nowext.out')
KTP_LOG_STR = pathtools.read_file(OUT_PATH, 'dmm_nowext.log')

WELL_ENE_DCT = {'DMM-R2': 0.006215045609626311, 'FakeW-CH3+CH3OCHO': None,
                'DMM-R1': 0.01466113323296463, 'FakeW-CH2O+CH3OCH2': None}

def test__well_energies():
    """ test mess_io._wellextend.well_energies
        also calls: _get_well_reactions, _max_temp_well_exists,
        mess_io.reader.well_thermal_energy
    """

    well_ene_dct = mess_io.well_energies(KTP_OUT_STR, KTP_LOG_STR, 0.1)
    for key, val in well_ene_dct.items():
        if val:
            assert numpy.isclose(val, WELL_ENE_DCT[key])
        else:
            assert WELL_ENE_DCT[key] is None
            # important to test this for fake wells

def test_format_well_extension_inp():
    CHECKLIST = ['WellExtension', 'ExtensionCorrection    0.6',
                 '  WellExtensionCap[kcal/mol]    3.90', '  WellExtensionCap[kcal/mol]    9.20']
    new_inp_str = mess_io._format_well_extension_inp(KTP_INP_STR, WELL_ENE_DCT, None)
    checklist = []
    for line in new_inp_str.splitlines():
        if 'Cap' in line or 'Extension' in line:
            checklist.append(line)
    
    assert CHECKLIST == checklist
    
""" missing tests: functions for well reduction (not really part of current workflow)
well_lumped_input_file (calls functions below, used in mechdriver)
well_lumping_scheme (generally None, so skips call to function below)
_format_well_extension_inp with well_lumping_scheme not None
"""

if __name__ == '__main__':
    test__well_energies()
    test_format_well_extension_inp()
