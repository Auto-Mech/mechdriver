"""
Sets up pivot point locations for a species
"""

import numpy as np
from qcelemental import constants as qcc

# Conversion factors
DEG2RAD = qcc.conversion_factor('degree', 'radian')
RAD2DEG = qcc.conversion_factor('radian', 'degree')


# def calc_pivot_point_internal(xyzp, xyz1):
#     """
#     """
#
#     dxyzp = xyzp - xyz1
#
#     # rho
#     rho = np.linalg.norm(dxyzp)
#
#     # theta
#     if abs(dxyzp[0]) < 0.01 and abs(dxyzp[1]) < 0.01:
#         theta = 0.01
#     elif abs(dxyzp[0]) < 0.01 and abs(dxyzp[1]) < 0.01:
#         theta = (np.pi / 2.0) * RAD2DEG
#     else:
#         theta = np.atan2(dxyzp[1], dxyzp[0]) * (180.0 / np.pi)
#
#     # phi
#     if abs(dxyzp[2]) < 0.01:
#         phi = (np.pi / 2.0) * RAD2DEG
#     else:
#         phi = ((atan2(np.sqrt(dxyzp[0]**2 + dxyzp[1]**2), dxyzp[2])) *
#                (180.0 / np.pi)
#
#     return rho, theta, phi



def calc_pivot_point_xyz(xyz1, xyz2, xyz3, in_dist, in_angle, in_dihed):
    """ geometric approach for calculating the xyz coordinates of atom A
        when the xyz coordinates of the A B and C are known and
        the position is defined w/r to A B C with internal coordinates
    """

    # Build initial coordinates
    #xyz1 = np.array(xyz1, dtype=float)
    #xyz2 = np.array(xyz2, dtype=float)
    #xyz3 = np.array(xyz3, dtype=float)
    in_angle *= DEG2RAD
    in_dihed *= DEG2RAD

    # Get the A, B, C, and P points in the rt system
    xyzp_rt = np.array([in_dist * np.sin(in_angle) * np.cos(in_dihed),
                        in_dist * np.cos(in_angle),
                        -(in_dist * np.sin(in_angle) * np.sin(in_dihed))])
    # xyz1_rt = np.array([0.0, 0.0, 0.0])

    dist12 = np.linalg.norm(xyz1 - xyz2)
    dist13 = np.linalg.norm(xyz1 - xyz3)
    dist23 = np.linalg.norm(xyz2 - xyz3)
    val = ((dist12**2 + dist13**2 - dist23**2) / 2.0 / dist12)
    xyz2_rt = np.array([0.0, dist12, 0.0])
    xyz3_rt = np.array([val, np.sqrt(dist13**2 - val**2), 0.0])

    # calculate the check1?
    # xyzp_rt_tmp = np.array([xyzp_rt[0], xyzp_rt[1], -xyzp_rt[2]])
    # check1 = np.linalg.norm(xyzp_rt_tmp - xyz3_rt)

    # translate original frame of ref coors so that xyz1 is at (0, 0, 0)
    # xyz1_t = np.array([0.0, 0.0, 0.0])
    xyz2_t = xyz2 - xyz1
    xyz3_t = xyz3 - xyz1

    # rotation matrix to rotate back to the original ref system
    r12 = (xyz2[0] - xyz1[0]) / xyz2_rt[1]
    r22 = (xyz2[1] - xyz1[1]) / xyz2_rt[1]
    r32 = (xyz2[2] - xyz1[2]) / xyz2_rt[1]

    r11 = (xyz3[0] - xyz1[0] - xyz3_rt[1]*r12) / xyz3_rt[0]
    r21 = (xyz3[1] - xyz1[1] - xyz3_rt[1]*r22) / xyz3_rt[0]
    r31 = (xyz3[2] - xyz1[2] - xyz3_rt[1]*r32) / xyz3_rt[0]

    anum_aconst = (xyz2_t[1] - xyz3_t[1]) / (xyz3_t[0]*xyz2_t[0])
    den_aconst = (xyz2_t[2] - xyz3_t[2]) / (xyz3_t[0]*xyz2_t[0])

    if abs(anum_aconst) < 1.0e-6 and abs(den_aconst) < 1.0e-6:
        if anum_aconst < 0.0:
            aconst = -1.0e20
        else:
            aconst = 1.0e20
    elif abs(den_aconst) < 1.0e-6:
        if anum_aconst < 0.0:
            aconst = -1.0e20
        else:
            aconst = 1.0e20
    else:
        aconst = (((xyz2_t[1] - xyz3_t[1]) / (xyz3_t[0]*xyz2_t[0])) /
                  ((xyz2_t[2] - xyz3_t[2]) / (xyz3_t[0]*xyz2_t[0])))

    den1 = (xyz3_t[1] / xyz3_t[0]) - aconst * (xyz3_t[2] / xyz3_t[0])
    if den1 == 0.0:
        den1 = 1.0e-20
    bconst = 1.0 / den1

    # Set vals for another point
    valx = -(1.0 / np.sqrt(1.0 + (bconst**2) * (1.0 + aconst**2)))
    valy = -(valx * bconst)
    xyz4_t = np.array([valx, valy, -(valy * aconst)])

    r13 = xyz4_t[0]
    r23 = xyz4_t[1]
    r33 = xyz4_t[2]

    r13n = -r13
    r23n = -r23
    r33n = -r33

    # now rotate and translate back
    # here I check  the (001) vector direction to decide whether
    # to take the positive of negative results of the
    # square root taken above
    xap = (xyz1[0] + (r11 * xyzp_rt[0]) +
           (r12 * xyzp_rt[1]) + (r13 * xyzp_rt[2]))
    yap = (xyz1[1] + (r21 * xyzp_rt[0]) +
           (r22 * xyzp_rt[1]) + (r33 * xyzp_rt[2]))
    zap = (xyz1[2] + (r31 * xyzp_rt[0]) +
           (r32 * xyzp_rt[1]) + (r33 * xyzp_rt[2]))

    xan = (xyz1[0] + (r11 * xyzp_rt[0]) +
           (r12 * xyzp_rt[1]) + (r13n * xyzp_rt[2]))
    yan = (xyz1[1] + (r21 * xyzp_rt[0]) +
           (r22 * xyzp_rt[1]) + (r23n * xyzp_rt[2]))
    zan = (xyz1[2] + (r31 * xyzp_rt[0]) +
           (r32 * xyzp_rt[1]) + (r33n * xyzp_rt[2]))

    bvec = xyz1 - xyz2
    cvec = xyz2 - xyz3
    # np.array([(xyz1[0] - xyz2[1]),
    #                 (xyz1[1] - xyz2[1]),
    #                 (xyz1[2] - xyz2[2])])
    # cvec = np.array([(xyz2[0] - xyz3[1]),
    #                 (xyz2[1] - xyz3[1]),
    #                 (xyz2[2] - xyz3[2])])

    vec1 = (bvec[1] * cvec[2]) - (bvec[2] * cvec[1])
    vec2 = (bvec[2] * cvec[0]) - (bvec[0] * cvec[2])
    vec3 = (bvec[0] * cvec[1]) - (bvec[1] * cvec[0])

    if abs(xyz4_t[0]) > 1.0e-5:
        checkv = vec1 / xyz4_t[0]
    elif abs(xyz4_t[1]) > 1.0e-5:
        checkv = vec2 / xyz4_t[1]
    else:
        checkv = vec3 / xyz4_t[2]

    if checkv >= 0.0:
        xyzp = np.array([xap, yap, zap])
    else:
        xyzp = np.array([xan, yan, zan])

    return xyzp[0], xyzp[1], xyzp[2]


if __name__ == '__main__':
    XYZ1 = np.array([-3.1610493755, -0.1014051234, 0.0250282314])
    XYZ2 = np.array([-3.4995568399, 0.3704026023, -0.9202655294])
    XYZ3 = np.array([-3.5632022224, 0.4687104654, 0.8875892033])
    IN_DIST = 1.109
    IN_ANGLE = 109.5
    IN_DIHED = -120.0
    XP, YP, ZP = calc_pivot_point_xyz(
        XYZ1, XYZ2, XYZ3, IN_DIST, IN_ANGLE, IN_DIHED)
    print('\ncalculated:')
    print(XP, YP, ZP)
    print('\nactual:')
    print(-2.0523345821, -0.0976810188, 0.0638579578)
