""" cartesian geometry writers
"""


def write(symbs, xyzs):
    """ Write a molecular geometry to a string:
           symb1  xyz1 xyz2 xyz3
           symbn  xyzn xyzn xyzn

        :param symbs: atomic symbols of the atoms
        :type symbs: tuple(str)
        :param xyzs: xyz coordinates of the atoms
        :type xyzs: tuple(float)
        :rtype: str
    """

    natms = len(symbs)
    assert len(xyzs) == natms

    geo_str = '\n'.join(
        f'{symb:2s} {xyz[0]:10.6f} {xyz[1]:10.6f} {xyz[2]:10.6f}'
        for symb, xyz in zip(symbs, xyzs))

    return geo_str


def write_xyz(symbs, xyzs, comment=None):
    """ Write a molecular geometry to a string:
           natom
           comment
           symb1  xyz1 xyz2 xyz3
           symbn  xyzn xyzn xyzn

        :param symbs: atomic symbols of the atoms
        :type symbs: tuple(str)
        :param xyzs: xyz coordinates of the atoms
        :type xyzs: tuple(float)
        :param comment: string to place in the comment line of string
        :type comment: str
        :rtype: str
    """

    comment = '' if comment is None else comment

    natms = len(symbs)
    assert len(xyzs) == natms

    geo_str = write(symbs=symbs, xyzs=xyzs)
    xyz_str = f' {natms:d}\n{comment:s}\n{geo_str:s}'

    return xyz_str


def write_xyz_trajectory(symbs, xyzs_lst, comments=None):
    """ Write a series of molecular geometries to trajectory file which
        is a string that collated by several xyz-file format geometry strings.

        :param symbs: atomic symbols of the atoms
        :type symbs: tuple(str)
        :param xyzs_lst: xyz coordinates for a set of molecular geometries
        :type xyzs_lst: tuple(tuple(float))
        :param comments: list of comments for each of the molecular geometries
        :type comments: tuple(str)
        :rtype: str
    """

    ngeos = len(xyzs_lst)
    comments = ('',) * ngeos if comments is None else comments
    assert len(comments) == ngeos
    xyz_strs = [write_xyz(symbs, xyzs, comment)
                for xyzs, comment in zip(xyzs_lst, comments)]

    return '\n'.join(xyz_strs)
