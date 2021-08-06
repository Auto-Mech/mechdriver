""" utilites
"""


def nsamp_init(nsamp_par, ntaudof):
    """ determine nsamp for given species"""
    if nsamp_par[0]:
        nsamp = min(nsamp_par[1] + nsamp_par[2] * nsamp_par[3]**ntaudof,
                    nsamp_par[4])
        # print('Setting nsamp using formula: min(A+B*C**n')
    else:
        nsamp = nsamp_par[5]
    return nsamp
