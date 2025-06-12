"""
 tests chemkin_io.writer.transport.properties
"""

import automol.chi
from chemkin_io.writer.transport import properties as writer


SPC_TRANS_DCT = {
    'He': {'geo': automol.chi.geometry('InChI=1S/He'),
           'epsilon': 0.0000362344204,
           'sigma': 5.1306064,
           'dipole_moment': 0.00,
           'polarizability': 9.2898,
           'zrot': 2.00
           },
    'O2': {'geo': automol.chi.geometry('InChI=1S/O2/c1-2'),
           'epsilon': 0.00048935,
           'sigma': 3.458,
           'dipole_moment': 0.00,
           'polarizability': 10.7974,
           'zrot': 3.800
           },
    'H2O': {'geo': automol.chi.geometry('InChI=1S/H2O/h1H2'),
            'epsilon': 0.002017423107,
            'sigma': 5.561464,
            'dipole_moment': 1.851,
            'polarizability': 64.07556,
            'zrot': 4.00
            }
}

TRANS_STR = """! THEORETICAL TRANSPORT PROPERTIES
!
! (1) Shape, index denotes atom (0), linear molec. (1), nonlinear molec. (2);
! (2) Epsilon, the Lennard-Jones well depth, in K;
! (3) Sigma, the Lennard-Jones collision diameter, in Angstrom;
! (4) Mu, total dipole moment, in Debye;
! (5) Alpha, mean static polarizability, in Angstrom^3; and
! (6) Z_rot, rotational relaxation collision number at 298 K.
! Species           Shape     Epsilon   Sigma      Mu   Alpha   Z_Rot
He                      0      11.442   2.715   0.000   1.377   2.000
O2                      1     154.525   1.830   0.000   1.600   3.800
H2O                     2     637.052   2.943   1.851   9.495   4.000
"""


def test__transport_write():
    """ Tests the transport writer
    """
    trans_str = writer(SPC_TRANS_DCT)
    assert trans_str == TRANS_STR
