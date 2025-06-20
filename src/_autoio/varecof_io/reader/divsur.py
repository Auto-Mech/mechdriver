"""
  Functions to parse information out of the divsur input and output files
"""

import automol.geom


def frame_geometries(divsur_str):
    """ Read the geometries of the fragments out of divsur.out file string
        which are in the coord system defined in the divsur.inp frames.

        :param divsur_str: string with geoms made from divsur frames
        :type divsur_str: str
        :return fgeo1: geometry of fragment 1 in divsur frame
        :rtype: str
        :return fgeo2: geometry of fragment 2 in divsur frame
        :rtype: str
    """

    # Get where the fragment geometries are defined
    lines = divsur_str.splitlines()
    for i, line in enumerate(lines):
        if 'Fragment 1:' in line:
            f1_idx = i+1
        if 'Fragment 2:' in line:
            f2_idx = i+1
        if 'SURFACE INDEX' in line:
            end_idx = i
            break

    if f1_idx is not None and f2_idx is not None:
        fgeo1_str = '\n'.join(lines[f1_idx: f2_idx-2])
        fgeo2_str = '\n'.join(lines[f2_idx: end_idx])
        fgeo1 = automol.geom.from_string(fgeo1_str)
        fgeo2 = automol.geom.from_string(fgeo2_str)
    else:
        fgeo1 = None
        fgeo2 = None

    return fgeo1, fgeo2
