"""
  NEW: Handle rotational data info
"""

from mechlib.amech_io import printer as ioprinter


def read_geom(pf_filesystems):
    """ Read a Cartesian geometry from the SAVE filesystem for a
        species or transition state.

        something about electronic structure mehtod in pf_filesys
    """

    # Get the harmonic filesys information
    [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['harm']

    # Read the filesys for the geometry
    if cnf_path:
        geom = cnf_fs[-1].file.geometry.read(min_cnf_locs)
        ioprinter.reading('geometry', cnf_path)
    else:
        ioprinter.info_message('No geometry found at path:', cnf_path)

    return geom


def read_rotational_values(pf_filesystems):
    """ Read the vibration-rotation matrix and centrifugal distortion
        constant matrix calculated via VPT2 from the SAVE filesystem
        for a species or transition state.
    """

    # Set up vpt2 level filesystem for rotational values
    [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['vpt2']

    # Read the filesys for rotational anharmonicity information
    if cnf_path:
        ioprinter.reading(
            'Vib-Rot Matrix and Centrifugal Dist. Consts from path:', cnf_path)
        vibrot_mat = cnf_fs[-1].file.vibro_rot_alpha_matrix.read(
            min_cnf_locs)
        cd_consts = cnf_fs[-1].file.quartic_centrifugal_dist_consts.read(
            min_cnf_locs)
    else:
        ioprinter.info_message(
            'No Vib-Rot Matrix and Centrifugal Dist. Consts from path:',
            cnf_path)

    return vibrot_mat, cd_consts
