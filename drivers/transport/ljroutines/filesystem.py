"""
Routines for interacting with the autofile filesystem
"""

import autofile
import moldr


def read_lj_from_save(target_save_prefix, target_info, theory_level):
    """ Reads existing Lennard-Jones Sigma and Epsilon from
        the Save Filesystem
    """

    sigma = None
    epsilon = None

    # Set new theory level for the filesystem
    theory_level = [theory_level[1],
                    theory_level[2],
                    moldr.util.orbital_restriction(
                        target_info, theory_level)]

    # Search save file system for LJ params
    tgt_save_fs = autofile.fs.species(target_save_prefix)
    if tgt_save_fs.leaf.exists(target_info):
        tgt_save_path = tgt_save_fs.leaf.path(target_info)
        etrans_save_fs = autofile.fs.energy_transfer(tgt_save_path)
        if etrans_save_fs.leaf.exists(theory_level):
            sigma = etrans_save_fs.leaf.file.lennard_jones_sigma.read(
                theory_level)
            epsilon = etrans_save_fs.leaf.file.lennard_jones_epsilon.read(
                theory_level)

    return sigma, epsilon


def write_lj_to_save(sigma, epsilon,
                     target_save_prefix, target_info, theory_level):
    """ Writes the computed Lennard-Jones Sigma and Epsilon
        to the Save FileSystem
    """

    # Set new theory level for the filesystem
    theory_level = [theory_level[1],
                    theory_level[2],
                    moldr.util.orbital_restriction(
                        target_info, theory_level)]

    # Search save file system for LJ params
    tgt_save_fs = autofile.fs.species(target_save_prefix)
    tgt_save_fs.leaf.create(target_info)
    tgt_save_path = tgt_save_fs.leaf.path(target_info)
    etrans_save_fs = autofile.fs.energy_transfer(tgt_save_path)
    etrans_save_fs.leaf.create(theory_level)
    sigma = etrans_save_fs.leaf.file.lennard_jones_sigma.write(
        sigma, theory_level)
    epsilon = etrans_save_fs.leaf.file.lennard_jones_epsilon.write(
        epsilon, theory_level)

    etrans_save_path = etrans_save_fs.leaf.path(theory_level)
    print('\nWriting Lennard-Jones parameters to Save FileSystem at\n')
    print(etrans_save_path)
