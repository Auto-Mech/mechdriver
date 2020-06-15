"""
  NEW: Handle rotational data info
"""

def read_geom(pf_filesystems):
    """ Read the geometry from the filesys
    """

    # Get the harmonic filesys information
    [cnf_fs, cnf_path, min_cnf_locs, save_path, _] = pf_filesystems['harm']

    # Read the filesys for the geometry
    if min_cnf_locs:
        geom = cnf_fs[-1].file.geometry.read(min_cnf_locs)
        print('Reading geometry from path:')
        print(cnf_path)
    else:
        print('No geometry found at path:')
        print(cnf_path)
    return geom


def read_rotational_values(pf_filesystems):
    """ Read the rotational info from filesys
    """

    # Set up vpt2 level filesystem for rotational values
    [cnf_fs, cnf_path, min_cnf_locs, save_path] = pf_filesystems['vpt2']

    # Read the filesys for rotational anharmonicity information
    if min_cnf_locs:
        print('Reading Vib-Rot Matrix and Centrifugal Dist. Consts from path:')
        print(cnf_path)
        vibrot_mat = cnf_fs[-1].file.vibro_rot_alpha_matrix.read(
            min_cnf_locs)
        cd_consts = cnf_fs[-1].file.quartic_centrifugal_dist_consts.read(
            min_cnf_locs)
    else:
        print('No Vib-Rot Matrix and Centrifugal Dist. Consts from path:')
        print(cnf_path)

    return vibrot_mat, cd_consts
