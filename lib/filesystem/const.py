"""
  Deal with constraints
"""

import automol


def build_constraint_dct(zma, tors_names):
    """ Build a dictionary of constraints
    """
    constraint_names = [name
                        for name_lst in tors_names
                        for name in name_lst]
    constraint_names.sort(key=lambda x: int(x.split('D')[1]))
    zma_vals = automol.zmatrix.values(zma)
    constraint_dct = dict(zip(
        constraint_names,
        (round(zma_vals[name], 2) for name in constraint_names)
    ))

    return constraint_dct
