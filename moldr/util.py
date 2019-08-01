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
#            'job_options': ['verytight'],
#            'job_options': ['verytight'],
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

    return script_str, opt_script_str, kwargs, opt_kwargs


def orbital_restriction(mul, restrict_open_shell=False):
    """ orbital restriction logical
    """
    if restrict_open_shell:
        orb_restr = True
    else:
        orb_restr = (mul == 1)
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


def reference_geometry(ich, chg, mul, method, basis, orb_restr, prefix,
                       geom_dct):
    """ obtain reference geometry
    if data for reference method exists use that
    then geometry dictionary takes precedence
    if nothing else from inchi
    """
    spc_fs = autofile.fs.species(prefix)
    spc_fs.leaf.create([ich, chg, mul])
    spc_path = spc_fs.leaf.path([ich, chg, mul])

    thy_fs = autofile.fs.theory(spc_path)
    thy_fs.leaf.create([method, basis, orb_restr])
    thy_path = thy_fs.leaf.path([method, basis, orb_restr])

    if thy_fs.leaf.file.geometry.exists([method, basis, orb_restr]):
        thy_path = thy_fs.leaf.path([method, basis, orb_restr])
        print('getting reference geometry from', thy_path)
        geo = thy_fs.leaf.file.geometry.read([method, basis, orb_restr])
    else:
        if ich in geom_dct:
            print('getting reference geometry from geom_dct')
            geo = geom_dct[ich]
        else:
            print('getting reference geometry from inchi')
            geo = automol.inchi.geometry(ich)
    return geo


def min_energy_conformer_locators(save_prefix):
    """ locators for minimum energy conformer """
    cnf_save_fs = autofile.fs.conformer(save_prefix)
    cnf_alocs_lst = cnf_save_fs.leaf.existing()
    if cnf_alocs_lst:
        cnf_enes = [cnf_save_fs.leaf.file.energy.read(alocs)
                    for alocs in cnf_alocs_lst]
        min_cnf_alocs = cnf_alocs_lst[cnf_enes.index(min(cnf_enes))]
    else:
        min_cnf_alocs = None
    return min_cnf_alocs


def reaction_energy(save_prefix, rxn_ich, rxn_chg, rxn_mul, method, basis,
                    restrict_open_shell):
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


def reagent_energies(save_prefix, rgt_ichs, rgt_chgs, rgt_muls, method, basis,
                     restrict_open_shell):
    """ reagent energies """
    enes = []
    for rgt_ich, rgt_chg, rgt_mul in zip(rgt_ichs, rgt_chgs, rgt_muls):
        spc_save_fs = autofile.fs.species(save_prefix)
        spc_save_path = spc_save_fs.leaf.path([rgt_ich, rgt_chg, rgt_mul])

        orb_restr = orbital_restriction(rgt_mul, restrict_open_shell)
        thy_save_fs = autofile.fs.theory(spc_save_path)
        thy_save_path = thy_save_fs.leaf.path([method, basis, orb_restr])

        min_cnf_alocs = min_energy_conformer_locators(thy_save_path)
        cnf_save_fs = autofile.fs.conformer(thy_save_path)
        ene = cnf_save_fs.leaf.file.energy.read(min_cnf_alocs)
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
