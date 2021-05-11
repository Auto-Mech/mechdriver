""" Runs automech
"""

import os
import subprocess
import pytest

REPO_PATH = os.getcwd()
DB_PATH = os.path.join(REPO_PATH, 'tests')
RUN_PATH = os.path.dirname(__file__)
INP_PATH = os.path.join(RUN_PATH, 'inp')
RUN_DAT_PATH = os.path.join(INP_PATH, 'run.dat')
TEMPLATE_PATH = os.path.join(RUN_PATH, 'inp_template')
RUN_TEMP_PATH = os.path.join(TEMPLATE_PATH, 'run.dat')
CMD_LINE = (
    'python -u bin/automech.py {0} >& {0}/run.log & disown %1'.format(RUN_PATH)
)


def test__init_geom():
    """ test automech.py init_geom
    """
    tsks = '\tspc  init_geom    runlvl=lvl_scf  inplvl=lvl_scf'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__conf_samp():
    """ test automech.py conf_samp
    """

    tsks = '\n'.join(['\tspc  init_geom    runlvl=lvl_scf  inplvl=lvl_scf',
                     '\tspc  conf_samp    runlvl=lvl_scf  inplvl=lvl_scf'])
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))
         
    subprocess.call(CMD_LINE.split())


def test__conf_opt():
    """ test automech.py conf_energy
    """

    tsks = '\n'.join(['\tspc  conf_opt    runlvl=lvl_scf  inplvl=lvl_scf  cnf_range=n2',
                      '\tspc  conf_opt    runlvl=lvl_scf  inplvl=lvl_scf  cnf_range=e2'])
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__conf_grad():
    """ test automech.py conf_grad
    """

    tsks = '\tspc  conf_grad    runlvl=lvl_scf  inplvl=lvl_scf'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__conf_hess():
    """ test automech.py conf_hess
    """

    tsks = '\tspc  conf_hess    runlvl=lvl_scf  inplvl=lvl_scf'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__conf_energy():
    """ test automech.py conf_energy
    """

    tsks = '\n'.join(['\tspc  conf_energy    runlvl=lvl_mp2  inplvl=lvl_scf  cnf_range=n2',
                      '\tspc  conf_energy    runlvl=lvl_mp2  inplvl=lvl_scf  cnf_range=e2'])
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__hr_scan():
    """ test automech.py hr_scan
    """

    tsks = '\tspc  hr_scan    runlvl=lvl_scf  inplvl=lvl_scf tors_model=1dhrfa'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__hr_grad():
    """ test automech.py hr_scan
    """

    tsks = '\tspc  hr_grad    runlvl=lvl_scf  inplvl=lvl_scf tors_model=1dhrfa'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__hr_hess():
    """ test automech.py hr_scan
    """

    tsks = '\tspc  hr_hess    runlvl=lvl_scf  inplvl=lvl_scf tors_model=1dhrfa'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__hr_energy():
    """ test automech.py hr_scan
    """

    tsks = '\tspc  hr_energy    runlvl=lvl_scf  inplvl=lvl_scf tors_model=1dhrfa'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__find_ts():
    """ test automech.py init_geom
    """
    tsks = '\tts  find_ts    runlvl=lvl_scf  inplvl=lvl_scf'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    with pytest.raises(TypeError):
        subprocess.call(CMD_LINE.split())


def test__ts_conf_samp():
    """ test automech.py conf_samp
    """

    tsks = '\tts  conf_samp    runlvl=lvl_scf  inplvl=lvl_scf'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))
         
    subprocess.call(CMD_LINE.split())


def test__ts_conf_grad():
    """ test automech.py conf_grad
    """

    tsks = '\tts  conf_grad    runlvl=lvl_scf  inplvl=lvl_scf'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__ts_conf_hess():
    """ test automech.py conf_hess
    """

    tsks = '\tts  conf_hess    runlvl=lvl_scf  inplvl=lvl_scf'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__ts_conf_energy():
    """ test automech.py conf_energy
    """

    tsks = '\tts  conf_energy    runlvl=lvl_mp2  inplvl=lvl_scf'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__ts_hr_scan():
    """ test automech.py hr_scan
    """

    tsks = '\tts  hr_scan    runlvl=lvl_scf  inplvl=lvl_scf tors_model=1dhrfa'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__ts_hr_grad():
    """ test automech.py hr_scan
    """

    tsks = '\tts  hr_grad    runlvl=lvl_scf  inplvl=lvl_scf tors_model=1dhrfa'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__ts_hr_hess():
    """ test automech.py hr_scan
    """

    tsks = '\tts  hr_hess    runlvl=lvl_scf  inplvl=lvl_scf tors_model=1dhrfa'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__ts_hr_energy():
    """ test automech.py hr_scan
    """

    tsks = '\tts hr_energy    runlvl=lvl_scf  inplvl=lvl_scf tors_model=1dhrfa'
    drivers = '\tes'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__write_messrate():
    """ test automech.py hr_scan
    """

    tsks = ''
    drivers = '\twrite_messrate'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__run_messrate():
    """ test automech.py hr_scan
    """

    tsks = ''
    drivers = '\trun_messrate'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


def test__run_fits():
    """ test automech.py hr_scan
    """

    tsks = ''
    drivers = '\trun_fits'
    # Format the run.dat with the run-save dir paths
    with open(RUN_TEMP_PATH, 'r') as file_obj:
        run_str = file_obj.read()
    with open(RUN_DAT_PATH, 'w') as file_obj:
        file_obj.write(run_str.format(DB_PATH, tsks, drivers))

    subprocess.call(CMD_LINE.split())


if __name__ == '__main__':
 #   test__init_geom()
 #   test__conf_samp()
 #   test__conf_opt()
 #   test__conf_grad()
 #   test__conf_hess()
 #   test__conf_energy()
 #   test__hr_scan()
 #   test__hr_grad()
 #   test__hr_hess()
 #   test__hr_energy()

    test__find_ts()
 #   test__ts_conf_samp()
 #   test__ts_conf_grad()
 #   test__ts_conf_hess()
 #   test__ts_conf_energy()
 #   test__ts_hr_scan()
 #   test__ts_hr_grad()
 #   test__ts_hr_hess()
 #   test__ts_hr_energy()

 #   test__write_messrate()
 #   test__run_messrate()
 #   test__run_fits()
  
