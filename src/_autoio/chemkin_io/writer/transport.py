"""
 Writes the string for a Chemkin transport file
"""

from phydat import phycon
from chemkin_io.writer import _util as util


def properties(spc_trans_dct):
    """ Writes the string in containing data from several mechanism species
        used in calculating transport properties during Chemkin simulations.

        :param spc_trans_dct:
        # below is output units, not input
        :type spc_trans_dict: {spc_name:
            {'geo: molecular geometry object (ang)
             'epsilon': Lennard-Jones well depth epsilon/k_B (Kelvin),
             'sigma': Lennard-Jones collision diameter (angstroms),
             'dipole_moment': dipole moment (Debye),
             'polarizability': polarizability (angstroms^3),
             'zrot': rotational relaxation collision number at 298 K
            }
        }
        :return: chemkin_str: Chemkin string with transport data
        :rtype: str
    """

    # Initialize string with common header
    chemkin_str = util.CKIN_TRANS_HEADER_STR
    chemkin_str += '\n'

    # Add the headers for each of the columns
    chemkin_str += (
        f'{"! Species":20s}' +
        f'{"Shape":>5s}' +
        f'{"Epsilon":>12s}' +
        f'{"Sigma":>8s}' +
        f'{"Mu":>8s}' +
        f'{"Alpha":>8s}' +
        f'{"Z_Rot":>8s}'
    )
    chemkin_str += '\n'

    # Add the values to the string, formatting as necessary
    for name, dct in spc_trans_dct.items():

        geo = dct.get('geo')
        shape_idx = util.format_shape_idx(geo) if geo is not None else 2

        eps = dct.get('epsilon', 0.00) * phycon.EH2K
        sig = dct.get('sigma', 0.00) * phycon.BOHR2ANG
        dmom = dct.get('dipole_moment', 0.00)
        polar = dct.get('polarizability', 0.00) * phycon.BOHR2ANG**3
        zrot = dct.get('zrot', 1.00)

        chemkin_str += (
            f'{name:20s}' +
            f'{shape_idx:>5d}' +
            f'{eps:>12.3f}' +
            f'{sig:>8.3f}' +
            f'{dmom:>8.3f}' +
            f'{polar:>8.3f}' +
            f'{zrot:>8.3f}'
        )
        chemkin_str += '\n'

    return chemkin_str
