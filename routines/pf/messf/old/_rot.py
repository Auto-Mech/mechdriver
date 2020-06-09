"""
  Handle rotational data info
"""


def read_rotational_values(vpt2_save_fs, vpt2_min_cnf_locs):
    """ Read the rotational info from filesys
    """
    if vpt2_min_cnf_locs is not None:
        vibrot_mat = vpt2_save_fs[-1].file.vibro_rot_alpha_matrix.read(
            vpt2_min_cnf_locs)
        cd_consts = vpt2_save_fs[-1].file.quartic_centrifugal_dist_consts.read(
            vpt2_min_cnf_locs)

    return vibrot_mat, cd_consts
