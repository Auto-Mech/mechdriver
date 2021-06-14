"""
  NEW: Handle flux files
"""

import autofile


def read_flux(ts_save_path, vrc_locs=(0,)):
    """ Read the geometry from the filesys
    """

    vrc_fs = autofile.fs.vrctst(ts_save_path)
    if vrc_fs[-1].file.flux.exists(vrc_locs):
        flux_str = vrc_fs[-1].file.flux.read(vrc_locs)
    else:
        flux_str = None

    return flux_str
