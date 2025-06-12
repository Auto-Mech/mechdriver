""" Functions that take data parsed from INTDER output
    as well as other functionality and calculate useful
    values
"""

from automol.zmat import distance_coordinate_name
from automol.zmat import central_angle_coordinate_name
from automol.zmat import dihedral_angle_coordinate_name


def ted_zmatrix_coordinates(zma, mode_idx, intl_coords, ted_assignments):
    """ Take the output of a Total Energy Distribution analysis
        from INTDER and determine what internal coordinates
        defined in a Z-Matrix make up that mode.

        The output has the internal coordinates defined in terms of
        atom indices of a Cartesian geometry+Hessian. Therefore, the
        input zmatrix must have an atom ordering that aligns with
        that Cartesian system used to generate the TED analysis.

        See intder_io.reader documentation for what intl_coords and
        ted_assignment variable objects are structured like.

        ret dct = {'R2': 68.9, 'A5': 11.4, 'D6': 10.0}

        :param zma: Z-Matrix with coordinates to get names for
        :typ zma: automol.zmat object
        :param intl_coords: definition of internal coordinates in INTDER
            in terms of the atom indices
        :type inl_coords: tuple(tuple(str, tuple(int)))
    """

    # Get the intl coord idxs associated with vibrational mode
    ted_mode_dct_in_ted_idxs = ted_assignments[mode_idx][1]

    # For each internal coordinate for the TED,
    # Get the type and atom and use this to get zmat coordinate name
    ted_mode_dct_in_zma_names = {}
    for intl_idx, intl_pct in ted_mode_dct_in_ted_idxs.items():
        typ, idxs = intl_coords[abs(intl_idx)]
        if typ == 'STRE':
            name = distance_coordinate_name(zma, *idxs)
        elif typ == 'BEND':
            name = central_angle_coordinate_name(zma, *idxs)
        elif typ == 'TORS':
            name = dihedral_angle_coordinate_name(zma, *idxs)
        ted_mode_dct_in_zma_names[name] = intl_pct

    return ted_mode_dct_in_zma_names


def ted_coordinate_indices(mode_idx, intl_coords, ted_assignments):
    """ Take the output of a Total Energy Distribution analysis
        from INTDER and write the indices of the internal coordinates
        that make the mode which were parsed from the INTDER output.

        The output has the internal coordinates defined in terms of
        atom indices of a Cartesian geometry+Hessian. Therefore, the
        input zmatrix must have an atom ordering that aligns with
        that Cartesian system used to generate the TED analysis.

        See intder_io.reader documentation for what intl_coords and
        ted_assignment variable objects are structured like.

        ret dct = {(1, 2): 68.9, '(5, 6, 7): 11.4, (8, 9, 10, 11): 10.0}

        :param zma: Z-Matrix with coordinates to get names for
        :typ zma: automol.zmat object
        :param intl_coords: definition of internal coordinates in INTDER
            in terms of the atom indices
        :type inl_coords: tuple(tuple(str, tuple(int)))
    """

    # Get the intl coord idxs associated with vibrational mode
    ted_mode_dct_in_ted_idxs = ted_assignments[mode_idx][1]

    # For each internal coordinate for the TED,
    # Get the type and atom and use this to get zmat coordinate name
    ted_mode_dct_in_xidxs = {}
    for intl_idx, intl_pct in ted_mode_dct_in_ted_idxs.items():
        _, idxs = intl_coords[abs(intl_idx)]
        ted_mode_dct_in_xidxs[frozenset(idxs)] = intl_pct

    return ted_mode_dct_in_xidxs
