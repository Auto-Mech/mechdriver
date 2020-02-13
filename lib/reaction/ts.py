"""
Find a TS from the grid as well as associated vdW wells
"""

import automol


def check_angle(ts_zma, dist_info, rxn_class):
    """ Check the angle to amend the dct
    """
    angle = None
    dist_name = dist_info[0]
    if 'abstraction' in rxn_class or 'addition' in rxn_class:
        brk_name = dist_info[3]
        if dist_name and brk_name:
            ts_bnd = automol.zmatrix.bond_idxs(
                ts_zma, dist_name)
            brk_bnd = automol.zmatrix.bond_idxs(
                ts_zma, brk_name)
            ang_atms = [0, 0, 0]
            cent_atm = list(set(brk_bnd) & set(ts_bnd))
            if cent_atm:
                ang_atms[1] = cent_atm[0]
                for idx in brk_bnd:
                    if idx != ang_atms[1]:
                        ang_atms[0] = idx
                for idx in ts_bnd:
                    if idx != ang_atms[1]:
                        ang_atms[2] = idx

            geom = automol.zmatrix.geometry(ts_zma)
            angle = automol.geom.central_angle(
                geom, *ang_atms)

    return angle
