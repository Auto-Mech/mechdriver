"""
  Handle the generation of necessary paths for various things
"""

import os
import random
import autofile


def projrot(run_path):
    """ Path to run projrot
    """
    # Set up the filesys
    rand_int = random.randint(0, 9999999)

    bld_locs = ['PROJROT', rand_int]
    bld_save_fs = autofile.fs.build(run_path)
    bld_save_fs[-1].create(bld_locs)
    projrot_path = bld_save_fs[-1].path(bld_locs)

    print('Run path for ProjRot:')
    print(projrot_path)

    return projrot_path


def mess_tors(run_path):
    """ Path to run MESS for torsion
    """

    # Set up the filesys
    bld_locs = ['PF', 0]
    bld_save_fs = autofile.fs.build(run_path)
    bld_save_fs[-1].create(bld_locs)
    pf_path = bld_save_fs[-1].path(bld_locs)

    print('Run path for MESSPF:')
    print(pf_path)

    return pf_path


# Set paths to MESS jobs
def messrate_path(prefix, pes_formula, sub_pes_idx):
    """ Build a simple mess path using the run prefix
    """
    pes_str = '{}_{}'.format(pes_formula, sub_pes_idx)
    return os.path.join(prefix, 'MESSRATE', pes_str)


def messpf_path(prefix, inchi):
    """ Build a simple mess path using the run prefix
    """
    spc_formula = automol.inchi.formula_string(inchi)
    ich_key = automol.inchi.inchi_key(inchi)
    return os.path.join(prefix, 'MESSPF', spc_formula, ich_key)


# Write MESS files
def write_mess_file(mess_inp_str, dat_str_dct, mess_path,
                    filename='mess.inp'):
    """ Write MESS file
        add overwrite keyword to handle overwrite MESS input
    """
    ioprinter.messpf('write_file', mess_path)

    # Write the MESS file
    if not os.path.exists(mess_path):
        os.makedirs(mess_path)
    with open(os.path.join(mess_path, filename), 'w') as mess_file:
        mess_file.write(mess_inp_str)

    # Write all of the data files needed
    if dat_str_dct:
        ioprinter.info_message('Writing the MESS data files...')
        for fname, fstring in dat_str_dct.items():
            # dat_path = os.path.join(mess_path, fname)
            if fstring:
                data_file_path = os.path.join(mess_path, fname)
                ioprinter.writing('data file', data_file_path)
                with open(data_file_path, 'w') as file_obj:
                    file_obj.write(fstring)
    # ioprinter.info_message(' - WARNING: File will be overwriten.')
    # ioprinter.info_message('No additional MESS input file will be written.')


def write_cwd_pf_file(mess_str, inchi, fname='pf.inp'):
    """ Write a copy of the MESS file in the current working directory
    """

    # Set starting path
    starting_path = os.getcwd()

    # Set the MESS paths and build dirs if needed
    jobdir_mess_path = messpf_path(starting_path, inchi)
    if not os.path.exists(jobdir_mess_path):
        os.makedirs(jobdir_mess_path)

    # Write the files
    file_path = os.path.join(jobdir_mess_path, fname)
    if os.path.exists(file_path):
        for i in range(1, 51):
            if not os.path.exists(file_path+str(i+1)):
                fin_path = file_path+str(i+1)
                break
    else:
        fin_path = file_path
    with open(fin_path, 'w') as file_obj:
        file_obj.write(mess_str)

    ioprinter.saving('MESS input copy', fin_path)

    return fin_path


def thermo_paths(spc_dct_i, run_prefix, idx):
    """ Set up the path for saving the pf input and output.
        Placed in a MESSPF, NASA dirs high in run filesys.
    """

    # Get the formula and inchi key
    spc_info = sinfo.from_dct(spc_dct_i)
    spc_formula = automol.inchi.formula_string(spc_info[0])
    ich_key = automol.inchi.inchi_key(spc_info[0])

    # PF
    bld_locs = ['PF', idx]
    bld_save_fs = autofile.fs.build(run_prefix)
    bld_save_fs[-1].create(bld_locs)
    bld_path = bld_save_fs[-1].path(bld_locs)
    ioprinter.debug_message(
        'preparing thermo for:', spc_info[0],
        bld_path, spc_formula, ich_key)
    spc_pf_path = os.path.join(bld_path, spc_formula, ich_key)

    # NASA polynomials
    bld_locs = ['NASA', idx]
    bld_save_fs = autofile.fs.build(run_prefix)
    bld_save_fs[-1].create(bld_locs)
    bld_path = bld_save_fs[-1].path(bld_locs)
    spc_nasa_path = os.path.join(bld_path, spc_formula, ich_key)

    return spc_pf_path, spc_nasa_path


def ckin_path():
    """ Set path and make directories
    """
    starting_path = os.getcwd()
    path = os.path.join(starting_path, 'ckin')
    if not os.path.exists(path):
        os.makedirs(path)

    return path
