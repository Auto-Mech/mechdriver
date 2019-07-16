""" utilites 
"""
import os
import stat
import subprocess
import warnings
import autofile

def species_theory_path(ich, chg, mult, method, basis, orb_restr, prefix):
    """ path to theory directory """
    spc_alocs = [ich, chg, mult]         # aloc = absolute locator
    thy_rlocs = [method, basis, orb_restr]  # rloc = relative locator
    thy_alocs = spc_alocs + thy_rlocs
    spc_afs = autofile.fs.species()
    thy_afs = autofile.fs.theory(spc_afs, 'species')
    thy_afs.theory.dir.create(prefix, thy_alocs)
    thy_path = thy_afs.theory.dir.path(prefix, thy_alocs)
    return thy_path


def reaction_theory_path(rxn_ichs, rxn_chgs, rxn_muls, ts_mul, method, basis, orb_restr, prefix):
    """ path to theory directory """
    rxn_alocs = [rxn_ichs, rxn_chgs, rxn_muls, ts_mul]
    thy_rlocs = [method, basis, orb_restr]
    thy_alocs = rxn_alocs + thy_rlocs
    rxn_afs = autofile.fs.reaction()
    thy_afs = autofile.fs.theory(rxn_afs, 'reaction')
    thy_afs.theory.dir.create(prefix, thy_alocs)
    thy_path = thy_afs.theory.dir.path(prefix, thy_alocs)
    return thy_path


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
        ene = cnf_afs.conf_trunk.file.energy.read(save_prefix, thy_alocs)
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

