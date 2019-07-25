""" utilites
"""
import os
import stat
import subprocess
import warnings
import autofile
import automol
import elstruct


def run_qchem_par(prog):
    """ dictionary of parameters for different electronic structure codes
    """

    if prog == 'g09':
        script_str = ("#!/usr/bin/env bash\n"
                      "g09 run.inp run.out >> stdout.log &> stderr.log")
        opt_script_str = script_str
        kwargs = {
            'memory': 20,
            'machine_options': ['%NProcShared=10'],
            'gen_lines': ['# int=ultrafine'],
        }
        opt_kwargs = {
            'memory': 20,
            'machine_options': ['%NProcShared=10'],
            'gen_lines': ['# int=ultrafine'],
            'feedback': True,
            # 'job_options': ['verytight'],
            'errors': [
                elstruct.Error.OPT_NOCONV
            ],
            'options_mat': [
                [{},
                 {},
                 {},
                 {'job_options': ['calcfc']},
                 {'job_options': ['calcfc']},
                 {'job_options': ['calcall']}]
            ],
        }

    if prog == 'psi4':
        script_str = ("#!/usr/bin/env bash\n"
                      "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")
        opt_script_str = script_str
        kwargs = {}
        opt_kwargs = {}

    if prog == 'molpro':
        script_str = ("#!/usr/bin/env bash\n"
                      "molpro -n 8 run.inp -o run.out >> stdout.log &> stderr.log")
        opt_script_str = ("#!/usr/bin/env bash\n"
                          "molpro --mppx -n 12 run.inp -o run.out >> stdout.log &> stderr.log")
        kwargs = {
            'memory': 50,
        }
        opt_kwargs = {
            'memory': 50,
        }

    if prog == 'qchem':
        script_str = ("#!/usr/bin/env bash\n"
                      "molpro -i run.inp -o run.out >> stdout.log &> stderr.log")
        opt_script_str = script_str
        kwargs = {
            'memory': 50,
        }
        opt_kwargs = {}

    if prog == 'cfour':
        script_str = ("#!/usr/bin/env bash\n"
                      "molpro -i run.inp -o run.out >> stdout.log &> stderr.log")
        opt_script_str = script_str
        kwargs = {
            'memory': 50,
        }
        opt_kwargs = {}

    if prog == 'orca':
        script_str = ("#!/usr/bin/env bash\n"
                      "molpro -i run.inp -o run.out >> stdout.log &> stderr.log")
        opt_script_str = script_str
        kwargs = {
            'memory': 50,
        }
        opt_kwargs = {}

    return  script_str, opt_script_str, kwargs, opt_kwargs

def orbital_restriction(mult, restrict_open_shell=False):
    """ orbital restriction logical
    """
    if restrict_open_shell:
        orb_restr = True
    else:
        orb_restr = (mult == 1)
    return orb_restr


def geometry_dictionary(geom_path):
    """ read in dictionary of saved geometries
    """
    geom_dct = {}
    for dir_path, _, file_names in os.walk(geom_path):
        for file_name in file_names:
            file_path = os.path.join(dir_path, file_name)
            if file_path.endswith('.xyz'):
                xyz_str = autofile.file.read_file(file_path)
                geo = automol.geom.from_xyz_string(xyz_str)
                ich = automol.geom.inchi(geo)
                if ich in geom_dct:
                    print('Warning: Dupilicate xyz geometry for ', ich)
                geom_dct[ich] = geo
    return geom_dct


def reference_geometry(ich, chg, mult, method, basis, orb_restr, prefix, geom_dct):
    """ obtain reference geometry
    if data for reference method exists use that
    then geometry dictionary takes precedence
    if nothing else from inchi
    """
    spc_path = species_path(ich, chg, mult, prefix)
    thy_afs = autofile.fs.theory()
    thy_alocs = [method, basis, orb_restr]
    if thy_afs.theory.file.geometry.exists(spc_path, thy_alocs):
        thy_path = thy_afs.theory.dir.path(spc_path, thy_alocs)
        print('getting reference geometry from', thy_path)
        geo = thy_afs.theory.file.geometry.read(spc_path, thy_alocs)
    else:
        if ich in geom_dct:
            print('getting reference geometry from geom_dct')
            geo = geom_dct[ich]
        else:
            print('getting reference geometry from inchi')
            geo = automol.inchi.geometry(ich)
    print(automol.geom.xyz_string(geo))
    return geo


def theory_path(method, basis, orb_restr, prefix):
    """ path to theory directory """
    thy_alocs = [method, basis, orb_restr]
    thy_afs = autofile.fs.theory()
    thy_afs.theory.dir.create(prefix, thy_alocs)
    thy_path = thy_afs.theory.dir.path(prefix, thy_alocs)
    return thy_path


def species_path(ich, chg, mult, prefix):
    """ path to species directory """
    spc_alocs = [ich, chg, mult]         # aloc = absolute locator
    spc_afs = autofile.fs.species()
    spc_afs.species.dir.create(prefix, spc_alocs)
    spc_path = spc_afs.species.dir.path(prefix, spc_alocs)
    return spc_path


def reaction_path(rxn_ichs, rxn_chgs, rxn_muls, ts_mul, prefix):
    """ path to reaction directory """
    rxn_alocs = [rxn_ichs, rxn_chgs, rxn_muls, ts_mul]
    rxn_afs = autofile.fs.reaction()
    rxn_afs.reaction.dir.create(prefix, rxn_alocs)
    rxn_path = rxn_afs.reaction.dir.path(prefix, rxn_alocs)
    return rxn_path


def min_energy_conformer_locators(save_prefix):
    """ locators for minimum energy conformer """
    cnf_afs = autofile.fs.conformer()
    cnf_alocs_lst = cnf_afs.conf.dir.existing(save_prefix)
    if cnf_alocs_lst:
        cnf_enes = [cnf_afs.conf.file.energy.read(save_prefix, alocs)
                    for alocs in cnf_alocs_lst]
        min_cnf_alocs = cnf_alocs_lst[cnf_enes.index(min(cnf_enes))]
    else:
        min_cnf_alocs = None
    return min_cnf_alocs


def reaction_energy(save_prefix, rxn_ich, rxn_chg, rxn_mul, method, basis, restrict_open_shell):
    """ reaction energy """
    rct_ichs, prd_ichs = rxn_ich
    rct_chgs, prd_chgs = rxn_chg
    rct_muls, prd_muls = rxn_mul
    rct_enes = reagent_energies(
        save_prefix, rct_ichs, rct_chgs, rct_muls, method, basis,
        restrict_open_shell)
    prd_enes = reagent_energies(
        save_prefix, prd_ichs, prd_chgs, prd_muls, method, basis,
        restrict_open_shell)
    return sum(prd_enes) - sum(rct_enes)


def reagent_energies(save_prefix, rgt_ichs, rgt_chgs, rgt_muls, method, basis, restrict_open_shell):
    """ reagent energies """
    enes = []
    spc_afs = autofile.fs.species()
    thy_afs = autofile.fs.theory(spc_afs, 'species')
    cnf_afs = autofile.fs.conformer(thy_afs, 'theory')

    for rgt_ich, rgt_chg, rgt_mul in zip(rgt_ichs, rgt_chgs, rgt_muls):
        print(rgt_ich, rgt_mul)
        if restrict_open_shell:
            orb_restr = True
        else:
            orb_restr = (rgt_mul == 1)
        spc_alocs = [rgt_ich, rgt_chg, rgt_mul]
        thy_rlocs = [method, basis, orb_restr]
        thy_alocs = spc_alocs + thy_rlocs
        thy_save_path = thy_afs.theory.dir.path(save_prefix, thy_alocs)
        min_cnf_alocs = min_energy_conformer_locators(thy_save_path)
        ene = cnf_afs.conf.file.energy.read(save_prefix, thy_alocs + min_cnf_alocs)
        enes.append(ene)
    return enes


def run_script(script_str, run_dir):
    """ run a program from a script
    """

    script_name = 'build.sh'
    with _EnterDirectory(run_dir):
        # write the submit script to the run directory
        with open(script_name, 'w') as script_obj:
            script_obj.write(script_str)

        # make the script executable
        os.chmod(script_name, mode=os.stat(script_name).st_mode | stat.S_IEXEC)

        # call the program
        try:
            subprocess.check_call('./{:s}'.format(script_name))
        except subprocess.CalledProcessError as err:
            # if the program failed, continue with a warning
            warnings.warn("run failed in {}".format(run_dir))


class _EnterDirectory():

    def __init__(self, directory):
        assert os.path.isdir(directory)
        self.directory = directory
        self.working_directory = os.getcwd()

    def __enter__(self):
        os.chdir(self.directory)

    def __exit__(self, _exc_type, _exc_value, _traceback):
        os.chdir(self.working_directory)
