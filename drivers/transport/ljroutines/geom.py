"""
Routines for getting the initial geometry for OneDMin
"""

import automol
import moldr
import py1dmin.interface


def get_geometry(spc_info, thry_lvl, save_prefix,
                 geom_dct={}, conf='sphere'):
    """ get the geometry
    """

    # Check if spc in dct, if so, read from file
    if spc in geom_dct:
        geom = read_xyz_file()
        geom = opt(geo, method, basis)
    else:
        assert conf in ('minimum', 'spherical')
        # Read the file system
        if conf == 'minimum':
            geom = _grab_min_energy_conformer(cnf_save_fs)
        elif conf == 'spherical':
            geom = _grab_spherical_conformer(cnf_save_fs)

        # Grab the minimum/spherical structure from set of confs from RDKit
        if geom is None:
            confs = automol.inchi.conformers(ich, nconfs=300)
            if conf == 'minimum':
                geo = confs[0]
            elif conf == 'spherical':
                geo = py1dmin.interface.util.roundify_geometry(confs)
            # Optimize the geometry
            geo = opt(geo, method, basis)

    # Format the geom into xyz strings
    geo_str = automol.geom.string(geo)

    return geo_str


def _grab_min_energy_conformer(cnf_save_fs):
    """ Reads the minimum-energy conformer from the save FileSystem
    """

    min_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
    if min_locs:
        min_conf = cnf_save_fs.leaf.file.geometry.read(min_locs)
    else:
        min_conf = None

    return min_conf


def _grab_spherical_conformer(cnf_save_fs):
    """ Reads the conformer from the Save FileSystem that is the most
        spherical.
    """
    geoms = _grab_all_conformers(cnf_save_fs)
    if geoms is not None:
        round_geom = py1dmin.interface.util.roundify_geometry(geoms)
    else:
        round_geom = None

    return round_geom


def _grab_all_conformers(cnf_save_fs):
    """ Reads all conformers from the Save FileSystem
    """

    cnf_locs_lst = cnf_save_fs.leaf.existing()
    if cnf_locs_lst:
        cnf_geoms = [cnf_save_fs.leaf.file.geometry.read(locs)
                     for locs in cnf_locs_lst]
    else:
        cnf_geoms = None
    return cnf_geoms
