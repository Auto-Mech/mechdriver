""" Tests the parsing of reaction strings
"""

import os
import numpy as np
import ioformat
from chemkin_io.parser.reaction import get_rxn_param_dct
from chemkin_io.parser.reaction import get_pes_dct
from chemkin_io.parser.reaction import get_rxn_osclass_dct
from chemkin_io.parser.reaction import get_rxn_cmt_dct, get_rxn_strs_dct

PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

REACTION_BLOCK_COMMENTS = "\
C4H7ORvE4fmAB0 = C4H7O4H74fm1                                             1.000E+00     0.000        0  ! pes.subpes.channel  1.1.1\n \
C4H7ORvE4fmAB0 = C4H7O-kSV4fm                                             1.000E+00     0.000        0  ! pes.subpes.channel  1.1.2\n \
C4H7ORvE4fmAB0 = C4H6O-RvErx51 + H-TcYTcY                                 1.000E+00     0.000        0  ! pes.subpes.channel  1.1.3\n \
C4H7ORvE4fmAA0 = C4H7O4H74fm0                                             1.000E+00     0.000        0  ! pes.subpes.channel  1.1.4\n \
C4H7ORvE4fmAA0 = C4H7O-kSV4fm                                             1.000E+00     0.000        0  ! pes.subpes.channel  1.1.5\n \
C4H7ORvE4fmAA0 = C4H6O-RvErx50 + H-TcYTcY                                 1.000E+00     0.000        0  ! pes.subpes.channel  1.1.6\n \
C4H7O4H74fm0 = C3H4OALAD-Wv9FbZ + CH3                                     1.000E+00     0.000        0  ! pes.subpes.channel  1.1.7\n \
C4H7O4H74fm0 = C2H4OALD-UPQWKw + C2H3ALK-S58hH1                           1.000E+00     0.000        0  ! pes.subpes.channel  1.1.8\n \
C4H7O-kSV4fm = C2H4OALD-UPQWKw + C2H3ALK-S58hH1                           1.000E+00     0.000        0  ! pes.subpes.channel  1.1.9\n \
C4H8ORvEsWvAA0 + OH = C4H7ORvE4fmAA0 + H2O                                1.000E+00     0.000        0  ! pes.subpes.channel  2.1.1\n \
C4H8ORvEsWvAB + OH = C4H7ORvE4fmAB0 + H2O                                 1.000E+00     0.000        0  ! pes.subpes.channel  2.2.2\n \
C4H8ORvEsWvAA0 + HO2-S580KW = C4H7ORvE4fmAA0 + H2O2-S58pAY                 1.000E+00     0.000        0  ! pes.subpes.channel  3.1.1\n \
C4H8ORvEsWvAA0 + CH3 = C4H7ORvE4fmAA0 + CH4                                1.000E+00     0.000        0  ! pes.subpes.channel  4.1.1\n \
C4H8ORvEsWvAA0 + CH3O2RO2-2LTcwB = C4H7ORvE4fmAA0 + CH4O2QOOH-2LTWKw       1.000E+00     0.000        0  ! pes.subpes.channel  5.1.1\n \
C4H8ORvEsWvAA0 + CH3O-S58cwB = C4H7ORvE4fmAA0 + CH4O-S58WKw                1.000E+00     0.000        0  ! pes.subpes.channel  6.1.1\n \
C4H8ORvEsWvAA0 + Cl = C4H7ORvE4fmAA0 + HCl                                 1.000E+00     0.000        0  ! pes.subpes.channel  7.1.1\n \
"
REACTION_BLOCK_OKCOMM = "\
C4H7ORvE4fmAB0 = C4H7O4H74fm1                                             1.000E+00     0.000        0  # pes.subpes.channel  1.1.1\n \
C4H7ORvE4fmAB0 = C4H7O-kSV4fm                                             1.000E+00     0.000        0  # pes.subpes.channel  1.1.2\n \
C4H7ORvE4fmAB0 = C4H6O-RvErx51 + H-TcYTcY                                 1.000E+00     0.000        0  # pes.subpes.channel  1.1.3\n \
C4H7ORvE4fmAA0 = C4H7O4H74fm0                                             1.000E+00     0.000        0  # pes.subpes.channel  1.1.4\n \
C4H7ORvE4fmAA0 = C4H7O-kSV4fm                                             1.000E+00     0.000        0  # pes.subpes.channel  1.1.5\n \
C4H7ORvE4fmAA0 = C4H6O-RvErx50 + H-TcYTcY                                 1.000E+00     0.000        0  # pes.subpes.channel  1.1.6\n \
C4H7O4H74fm0 = C3H4OALAD-Wv9FbZ + CH3                                     1.000E+00     0.000        0  # pes.subpes.channel  1.1.7\n \
C4H7O4H74fm0 = C2H4OALD-UPQWKw + C2H3ALK-S58hH1                           1.000E+00     0.000        0  # pes.subpes.channel  1.1.8\n \
C4H7O-kSV4fm = C2H4OALD-UPQWKw + C2H3ALK-S58hH1                           1.000E+00     0.000        0  # pes.subpes.channel  1.1.9\n \
C4H8ORvEsWvAA0 + OH = C4H7ORvE4fmAA0 + H2O                                1.000E+00     0.000        0  # pes.subpes.channel  2.1.1\n \
C4H8ORvEsWvAB + OH = C4H7ORvE4fmAB0 + H2O                                 1.000E+00     0.000        0  # pes.subpes.channel  2.2.2\n \
C4H8ORvEsWvAA0 + HO2-S580KW = C4H7ORvE4fmAA0 + H2O2-S58pAY                 1.000E+00     0.000        0  # pes.subpes.channel  3.1.1\n \
C4H8ORvEsWvAA0 + CH3 = C4H7ORvE4fmAA0 + CH4                                1.000E+00     0.000        0  # pes.subpes.channel  4.1.1\n \
C4H8ORvEsWvAA0 + CH3O2RO2-2LTcwB = C4H7ORvE4fmAA0 + CH4O2QOOH-2LTWKw       1.000E+00     0.000        0  # pes.subpes.channel  5.1.1\n \
C4H8ORvEsWvAA0 + CH3O-S58cwB = C4H7ORvE4fmAA0 + CH4O-S58WKw                1.000E+00     0.000        0  # pes.subpes.channel  6.1.1\n \
C4H8ORvEsWvAA0 + Cl = C4H7ORvE4fmAA0 + HCl                                 1.000E+00     0.000        0  # pes.subpes.channel  7.1.1\n \
"
REACTION_BLOCK_PESUBPES_WDUP = "REACTIONS     CAL/MOLE     MOLES\n\
\n\
NC3H7 = IC3H7                    1.000E+00     1.000        1  # pes.subpes.channel.submech_prompt  1.1.1.RAD_DECO_NC3H\n\
\n\
C3H6 + H = NC3H7                 1.040E+49   -11.500    15359  # pes.subpes.channel.submech_prompt  1.1.2.RAD_DECO_NC3H\n\
  PLOG /             1.300E-03   7.990E+81   -23.161    22239 /\n\
  PLOG /             1.300E-03   1.850E+26    -5.830     3866 /\n\
  PLOG /             4.000E-02   4.240E+68   -18.427    19665 /\n\
  PLOG /             4.000E-02   2.820E+30    -6.490     5471 /\n\
  PLOG /             1.000E+00   1.040E+49   -11.500    15359 /\n\
  PLOG /             1.000E+00   3.780E+28    -5.570     5625 /\n\
  PLOG /             1.000E+01   6.200E+41    -8.892    14637 /\n\
  PLOG /             1.000E+01   1.460E+25    -4.280     5248 /\n\
  PLOG /             1.000E+02   4.220E+27    -4.390     9346 /\n\
  PLOG /             1.000E+02   1.000E-10     0.000        0 /\n\
\n\
C3H6 + H = IC3H7                 3.260E+61   -14.940    20161  # pes.subpes.channel.submech_prompt  1.1.3.RAD_DECO_IC3H\n\
  PLOG /             1.300E-03   1.350E+44   -10.680     8196 /\n\
  PLOG /             1.300E-03   2.170E+130  -32.580   136140 /\n\
  PLOG /             4.000E-02   2.110E+57   -14.230    15147 /\n\
  PLOG /             4.000E-02   2.250E+29    -5.840     4242 /\n\
  PLOG /             1.000E+00   3.260E+61   -14.940    20161 /\n\
  PLOG /             1.000E+00   1.060E+30    -5.630     5613 /\n\
  PLOG /             1.000E+01   5.300E+56   -13.120    20667 /\n\
  PLOG /             1.000E+01   6.110E+26    -4.440     5182 /\n\
  PLOG /             1.000E+02   1.110E+50   -10.800    20202 /\n\
  PLOG /             1.000E+02   2.730E+23    -3.260     4597 /\n\
\n\
C3H6 + H = C2H4 + CH3            2.670E+12     0.470     5431  # pes.subpes.channel.submech_prompt  1.1.4.\n\
  PLOG /             1.300E-03   1.540E+09     1.350     2542 /\n\
  PLOG /             1.300E-03   1.000E-10     0.000        0 /\n\
  PLOG /             4.000E-02   7.880E+10     0.870     3600 /\n\
  PLOG /             4.000E-02   1.000E-10     0.000        0 /\n\
  PLOG /             1.000E+00   2.670E+12     0.470     5431 /\n\
  PLOG /             1.000E+00   1.000E-10     0.000        0 /\n\
  PLOG /             1.000E+01   9.250E+22    -2.600    12898 /\n\
  PLOG /             1.000E+01   1.240E+05     2.520     3679 /\n\
  PLOG /             1.000E+02   1.320E+23    -2.420    16500 /\n\
  PLOG /             1.000E+02   2.510E+03     2.910     3981 /\n\
\n\
C2H4 + CH3 = NC3H7               7.670E+47   -11.170    22366  # pes.subpes.channel.submech_prompt  1.1.5.RAD_DECO_NC3H\n\
  PLOG /             1.300E-03   8.670E+48   -12.540    18206 /\n\
  PLOG /             1.300E-03   1.120E+43   -11.300    13080 /\n\
  PLOG /             4.000E-02   1.060E+49   -12.040    20001 /\n\
  PLOG /             4.000E-02   7.280E+39    -9.880    13164 /\n\
  PLOG /             1.000E+00   7.670E+47   -11.170    22366 /\n\
  PLOG /             1.000E+00   2.600E+33    -7.460    12416 /\n\
  PLOG /             1.000E+01   1.810E+45   -10.030    23769 /\n\
  PLOG /             1.000E+01   3.850E+27    -5.380    11455 /\n\
  PLOG /             1.000E+02   2.040E+40    -8.250    24214 /\n\
  PLOG /             1.000E+02   1.660E+21    -3.170    10241 /\n\
\n\
C3H8 + H = IC3H7 + H2            3.060E+04     2.800     3473  # pes.subpes.channel.submech_prompt  2.1.1.RAD_GEN_IC3H7\n\
  DUP\n\
C3H8 + H = IC3H7 + H2            5.730E+13     0.100    10600\n\
  DUP\n\
\n\
C3H8 + H = NC3H7 + H2            1.290E+04     2.930     5246  # pes.subpes.channel.submech_prompt  2.2.2.RAD_GEN_NC3H7\n\
  DUP\n\
C3H8 + H = NC3H7 + H2            1.470E+10     1.310    10300\n\
  DUP\n\
\n\
C3H8 + OH = IC3H7 + H2O          3.738E+05     2.305     -561  # pes.subpes.channel.submech_prompt  3.1.1.RAD_GEN_IC3H7\n\
\n\
C3H8 + OH = NC3H7 + H2O          6.865E+06     2.000      677  # pes.subpes.channel.submech_prompt  3.2.2.RAD_GEN_NC3H7\n\
END\n\
"
REACTION_BLOCK_OSCLASS = "REACTIONS\n\
!#[REACTIONCLASS][A1-R][ROPEN]\n \
C6H5=LC6H4+H	1.64E+16	0	88395	! LPM 02-05-2024 FIT FROM MADDEN 1997 JPCA 10.1021/JP970723D\n \
PLOG/ 0.5	1.64E+16	0	88395	/\n \
PLOG/ 10	7.50E+16	0	90868	/\n \
PLOG/ 300	1.53E+17	0	92071	/\n \
!#[ENDREACTIONCLASS][A1-R][ROPEN]\n \
\n \
!#[REACTIONCLASS][C4.DT-R][ADD_C2.T-M]\n \
C4H3+C2H2=>LC6H4+H	1.99E+15	0.00	13493 ! LPM 02-05-2024 FIT FROM MADDEN 1997 JPCA 10.1021/JP970723D	\n \
PLOG/ 0.5	1.99E+15	0.00	13493	/\n \
PLOG/ 10	3.76E+15	0.00	18929	/\n \
PLOG/ 300	3.01E+07	0.00	22759	/\n \
C4H3+C2H2=>CYC6H4+H	2.70E+35	-6.77	24304	 ! LPM 02-05-2024 FIT FROM MADDEN 1997 JPCA 10.1021/JP970723D\n \
PLOG/ 0.5	2.70E+35	-6.77	24304	/\n \
PLOG/ 10	2.96E+18	-1.97	17957	/\n \
PLOG/ 300	3.34E-12	4.14	6779	/\n \
!#[ENDREACTIONCLASS][C4.DT-R][ADD_C2.T-M]\n \
END\n\
"
NONINT_RXN_PRODUCTS_AND_OTHER_EXCEPTIONS = "REACTIONS\n\
C5H5CH3+C5H4CH3=>0.5C10H7CH3+0.5CH3+0.5H2 +0.5C12H8+1.5H2+0.5H    .3000E+13    .000  23000.0     \n\
CH3C6H4+C6H5C2H2=>.5C16H10+.5C14H10+H2+H+H   +5.00000E+012 +0.00000E+000 +0.00000E+000 \n\
A+B=>1.5C+H 1.0 1.0 1.0\n\
TIC4H7Q2-I=>ET12IC4-1OOH+OH                       +4.45000000E+009 +8.60000000E-001 +1.08000000E+004 \n\
CH3CHO+C2H5CO3<=>CH3CO+C2H5CO3H                   +1.91580000E+000 +3.64260000E+000 +5.64190000E+003 \n\
DUP \n\
CH3CHO+C2H5CO3<=>CH3CO+C2H5CO3H                   +9.81000000E-004 +4.32010000E+000 +2.63610000E+003 \n\
DUP\n\
    H+OH(+M)<=>H2O(+M)                                +2.50000000E+013 +2.34000000E-001 -1.14000000E+002 \n\
LOW /                                             +4.53000000E+021 -1.81000000E+000 +4.92000000E+002 /\n\
TROE /                           +7.30000000E-001 +1.00000000E-030 +1.00000000E+030 +1.00000000E+030 /\n\
AR / +0.4000 / HE / +0.5700 / O2 / +1.0000 / N2 / +1.0000 / H2 / +1.5000 / H2O / +13.0000 / \n\
O+O+M<=>O2+M                                      +6.16500000E+015 -5.00000000E-001 +0.00000000E+000 \n\
HE / +0.8300 / AR / +0.8300 / CO / +1.9000 / CH4 / +2.0000 / H2 / +2.5000 / C2H6 / +3.0000 / CO2 / +3.8000 / H2O / +12.0000 / \n\
H+O+M<=>OHV+M                                     +4.43000000E+014 +0.00000000E+000 +1.00000000E+004 \n\
AR / +0.4000 / O2 / +0.8000 / CO2 / +2.0000 / H2 / +2.0000 / CH4 / +2.0000 / H2O / +12.0000 / \n\
H2O2(+M)<=>OH+OH(+M)                              +2.00000000E+012 +9.00000000E-001 +4.87490000E+004 \n\
LOW /                                             +2.35301810E+024 -2.29293580E+000 +4.87434050E+004 /\n\
TROE /                                            +4.30000000E-001 +1.00000000E-030 +1.00000000E+030 /\n\
HE / +0.4400 / O2 / +0.7900 / CO2 / +1.0600 / N2 / +1.5000 / CO / +2.8000 / H2 / +3.7000 / H2O / +5.1000 / H2O2 / +5.2000 / \n\
TIC4H7Q2-I=>ET12IC4-1OOH+OH                       +4.45000000E+009 +8.60000000E-001 +1.08000000E+004 \n\
END\n"


def test_cmts():
    """ Test inline comments reading - random checks
    """
    keys_to_cmts = [
        (('C4H2', 'C6H5'), ('C10H7',), (None,)), #dup
        (('C6H5', 'C6H4C2H'), ('C14H10',), (None,)), #simple
        (('CH3C6H4', 'C6H5C2H5'), ('C6H5C2H4C6H5', 'CH3'), (None,)), # written as irrev, but with no bw
    ]
    cmts_to_check = [
        ' LPM 01-21-2021 UPDATE C3 ! DUPLICATE CAUSE THERE WERE 2 DIFFERENT RADICALS OF C10H7; HPLIM IS OLD CRECK !lumped C10H8tyl as CRECK  ',
        '', 
        '',]
        
    cmts_to_check = dict(zip(keys_to_cmts, cmts_to_check))
    # avoids calling reaction block (from parser.mechanism), so no line clean up is done
    ckin_str = ioformat.pathtools.read_file(DAT_PATH, 'pah_block.dat')
    rxn_cmts_dct = get_rxn_cmt_dct(ckin_str)
    # check comments
    for rxn, cmt in cmts_to_check.items():
        assert cmt == rxn_cmts_dct[rxn]['inline']

def test_rxn_strs():
    """ Test reaction string dictionary from file
    """       
    check_one = ['C4H2+C6H5=C10H7                               5.000E+12    0.000     5000.00 ! LPM 01-21-2021 UPDATE C3 ! DUPLICATE CAUSE THERE WERE 2 DIFFERENT RADICALS OF C10H7; HPLIM IS OLD CRECK !lumped C10H8tyl as CRECK  \nPLOG / 0.01\t1.25E+88\t-21.83\t48863/\t!800-1700\nPLOG / 0.1\t1.54E+81\t-19.48\t50784/\t!800-1800\nPLOG / 1\t1.28E+76\t-17.68\t55146/\t!800-2000\nPLOG / 10\t2.22E+58\t-12.44\t49974/\t!800-2500\nPLOG / 100\t1.02E+35\t-5.809\t39347/\t!800-2500\nDUPLICATE', 
                 'C4H2+C6H5=C10H7                      5.00E+12    0.000     5000.00 !2-naphthyl\nPLOG / 0.01\t1.83E+87\t-21.56\t48728/\nPLOG / 0.1\t4.57E+85\t-20.72\t54421/\nPLOG / 1\t7.82E+72\t-16.75\t53604/\nPLOG / 10\t3.74E+63\t-13.83\t55912/\nPLOG / 100\t2.55E+45\t-8.515\t52892/\nDUPLICATE\n\n!#[REACTIONCLASS][A1-R][ADD_C4.TT-M]\n!#[ENDREACTIONCLASS][A1-R][ADD_C4.TT-M]']
    ckin_str = ioformat.pathtools.read_file(DAT_PATH, 'pah_block.dat')
    rxn_strs_dct = get_rxn_strs_dct(ckin_str)
    
    assert rxn_strs_dct[(('C4H2', 'C6H5'), ('C10H7',), (None,))] == check_one
    
def test_rxn_osclass():
    """ Test reading opensmoke class
    """
    dct_class = {(('C6H5',), ('LC6H4', 'H'), (None,)): {'speciestype': 'A1-R', 'reactiontype': 'ROPEN'}, 
                 (('C4H3', 'C2H2'), ('LC6H4', 'H'), (None,)): {'speciestype': 'C4.DT-R', 'reactiontype': 'ADD_C2.T-M'}, 
                 (('C4H3', 'C2H2'), ('CYC6H4', 'H'), (None,)): {'speciestype': 'C4.DT-R', 'reactiontype': 'ADD_C2.T-M'}}
    
    os_class = get_rxn_osclass_dct(REACTION_BLOCK_OSCLASS)
    assert os_class == dct_class

def test_nonint_rxn_products():
    """ Test PLOG with comments
    """
    NAMES = {(('CH3C6H4', 'C6H5C2H2'), ('.5C16H10', '.5C14H10', 'H2', 'H', 'H'), (None,)),
        (('C5H5CH3', 'C5H4CH3'), ('0.5C10H7CH3', '0.5CH3', '0.5H2', '0.5C12H8', '1.5H2', '0.5H'), (None,)),
        (('A', 'B'), ('1.5C', 'H'), (None,)),
        (('TIC4H7Q2-I',), ('ET12IC4-1OOH', 'OH',), (None,)),
        (('CH3CHO', 'C2H5CO3'), ('CH3CO', 'C2H5CO3H',), (None,)),
        (('H2O2',), ('OH', 'OH'), ('(+M)',)),
        (('H', 'OH'), ('H2O',), ('(+M)',)),
        (('O', 'O'), ('O2',), ('+M',)),
        (('H', 'O'), ('OHV',), ('+M',)),
        (('TIC4H7Q2-I',), ('ET12IC4-1OOH', 'OH'), (None,)),
        }

    rxn_param_dct1 = get_rxn_param_dct(NONINT_RXN_PRODUCTS_AND_OTHER_EXCEPTIONS, 'cal/mole', 'moles')
    # 3 reactions found
    assert all(name in rxn_param_dct1.keys() for name in NAMES)
    assert all(name in NAMES for name in rxn_param_dct1.keys())
    
            
def test_arr():
    """ Tests the Arrhenius reader
    """

    ckin_str1 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2+M=OH+O+M     1.000E+15     0.000    25000\n'
        '     N2/1.400/   AR/1.000/   \n\n\n'  # note three spaces before \n
        'END\n\n')
    ckin_str2 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2=OH+O     1.000E+15     0.000    25000\n'
        'DUP\n'
        'H+O2=OH+O     1.000E+15     0.000    25000\n'
        'DUP\n\n\n'
        'END\n\n')
    ckin_str3 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2=OH+O     1.000E+15     0.000    25000\n'
        'DUP\n'
        'H+O2=OH+O     1.000E+15     0.000    27000\n'
        'DUP\n\n\n'
        'END\n\n')

    rxn_param_dct1 = get_rxn_param_dct(ckin_str1, 'cal/mole', 'moles')
    rxn_param_dct2 = get_rxn_param_dct(ckin_str2, 'cal/mole', 'moles')
    rxn_param_dct3 = get_rxn_param_dct(ckin_str3, 'cal/mole', 'moles')

    for params in rxn_param_dct1.values():  # should only be one rxn
        arr_tuples = params.arr
        arr_collid = params.arr_collid
        assert len(arr_tuples) == 1
        for arr_tuple in arr_tuples:
            assert np.allclose(arr_tuple, [1e15, 0, 25000])
        for spc, eff in arr_collid.items():
            assert spc in ('N2', 'AR')
            assert np.isclose(1.4, eff) or np.isclose(1.0, eff)

    for params in rxn_param_dct2.values():
        arr_tuples = params.arr
        assert len(arr_tuples) == 1  # removes duplicate bc EA is the same
        for arr_tuple in arr_tuples:
            assert np.allclose(arr_tuple, [2e15, 0, 25000])

    for params in rxn_param_dct3.values():
        arr_tuples = params.arr
        assert len(arr_tuples) == 2  # should be dup 
        assert np.allclose(arr_tuples[0], [1e15, 0, 25000])
        assert np.allclose(arr_tuples[1], [1e15, 0, 27000])

def test_plog2():
    """ Test PLOG with comments
    """
    rxn_param_dct1 = get_rxn_param_dct(REACTION_BLOCK_OSCLASS, 'cal/mole', 'moles')
    # 3 reactions found
    assert len(rxn_param_dct1.keys()) == 3
    for params in rxn_param_dct1.values():
        plog_dct = params.plog
        assert len(plog_dct.keys()) == 3 #3 sets of pressure
    params = rxn_param_dct1[(('C6H5',), ('LC6H4', 'H'), (None,))]
    plog_dct = params.plog
    for arr_tuples in plog_dct.values():
        for arr_tuple in arr_tuples:
            assert np.isclose(arr_tuple[1], 0.)
    
def test_plog():
    """ Tests the PLOG reader and writer
    """

    ckin_str1 = 'REACTIONS     CAL/MOLE     MOLES\n\n' \
        'H+O2=OH+O     1.000E+15     0.000    25000   ! Duplicates exist at ' \
        '1 atm (see below); only single 1-atm fit is written\n' \
        '    PLOG /1.000E-01   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+00   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+00   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+02   1.000E+15     0.000    25000 /\n\n\nEND\n\n'

    # This one is for checking the incorrect way of writing duplicate PLOGS
    ckin_str2 = 'REACTIONS     CAL/MOLE     MOLES\n\n' \
        'H+O2=OH+O     1.000E+15     0.000    25000   ! Duplicates exist at ' \
        '1 atm (see below); only single 1-atm fit is written\n' \
        '    PLOG /1.000E-01   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+00   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+00   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+02   1.000E+15     0.000    25000 /\nDUP\n' \
        'H+O2=OH+O     1.000E+15     0.000    25000   ! Duplicates exist at ' \
        '1 atm (see below); only single 1-atm fit is written\n' \
        '    PLOG /1.000E-01   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+00   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+00   1.000E+15     0.000    25000 /\n' \
        '    PLOG /1.000E+02   1.000E+15     0.000    25000 /\nDUP\n\n\nEND' \
        '\n\n'

    rxn_param_dct1 = get_rxn_param_dct(ckin_str1, 'cal/mole', 'moles')
    for params in rxn_param_dct1.values():
        plog_dct = params.plog
        for pressure, arr_tuples in plog_dct.items():
            for arr_tuple in arr_tuples:
                assert np.allclose(arr_tuple, [1e15, 0, 25000])
            if pressure == 1:  # this pressure should have duplicate Arrhenius
                assert len(arr_tuples) == 2

    # Check the duplicate case
    rxn_param_dct2 = get_rxn_param_dct(ckin_str2, 'cal/mole', 'moles')
    for params in rxn_param_dct2.values():
        plog_dct = params.plog
        plog_dups = params.plog_dups
        # Check the plog_dct
        for pressure, arr_tuples in plog_dct.items():
            for arr_tuple in arr_tuples:
                assert np.allclose(arr_tuple, [1e15, 0, 25000])
            if pressure == 1:  # this pressure should have duplicate Arrhenius
                assert len(arr_tuples) == 2
        # Check the plog_dups
        assert len(plog_dups) == 1
        dup_dct = plog_dups[0]
        for pressure, arr_tuples in dup_dct.items():
            for arr_tuple in arr_tuples:
                assert np.allclose(arr_tuple, [1e15, 0, 25000])
            if pressure == 1:  # this pressure should have duplicate Arrhenius
                assert len(arr_tuples) == 2


def test_cheb():
    """ Tests the Chebyshev reader and writer
    """

    ckin_str1 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+00     0.000        0\n'
        '    TCHEB/      500.00     2000.00 /\n'
        '    PCHEB/        0.03      100.00 /\n'
        '    CHEB /           6           4 /\n'
        '    CHEB /  -1.620E+01  -1.183E-01  -5.423E-02  -1.476E-02 /\n'
        '    CHEB /   2.578E+00   1.614E-01   7.359E-02   1.880E-02 /\n'
        '    CHEB /   1.068E-01  -7.235E-02  -2.733E-02  -3.778E-03 /\n'
        '    CHEB /   3.955E-02   1.207E-02   3.402E-04  -2.695E-03 /\n'
        '    CHEB /   8.557E-03   4.345E-03   3.670E-03   1.608E-03 /\n'
        '    CHEB /   8.599E-04  -1.758E-03  -7.502E-04   7.396E-07 /\n\n\n'
        'END\n\n')

    # For duplicates
    ckin_str2 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+00     0.000        0\n'
        '    TCHEB/      500.00     2000.00 /\n'
        '    PCHEB/        0.03      100.00 /\n'
        '    CHEB /           6           4 /\n'
        '    CHEB /  -1.620E+01  -1.183E-01  -5.423E-02  -1.476E-02 /\n'
        '    CHEB /   2.578E+00   1.614E-01   7.359E-02   1.880E-02 /\n'
        '    CHEB /   1.068E-01  -7.235E-02  -2.733E-02  -3.778E-03 /\n'
        '    CHEB /   3.955E-02   1.207E-02   3.402E-04  -2.695E-03 /\n'
        '    CHEB /   8.557E-03   4.345E-03   3.670E-03   1.608E-03 /\n'
        '    CHEB /   8.599E-04  -1.758E-03  -7.502E-04   7.396E-07 /\n'
        'DUP\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+00     0.000        0\n'
        '    TCHEB/      500.00     2000.00 /\n'
        '    PCHEB/        0.03      100.00 /\n'
        '    CHEB /           6           4 /\n'
        '    CHEB /  -1.620E+01  -1.183E-01  -5.423E-02  -1.476E-02 /\n'
        '    CHEB /   2.578E+00   1.614E-01   7.359E-02   1.880E-02 /\n'
        '    CHEB /   1.068E-01  -7.235E-02  -2.733E-02  -3.778E-03 /\n'
        '    CHEB /   3.955E-02   1.207E-02   3.402E-04  -2.695E-03 /\n'
        '    CHEB /   8.557E-03   4.345E-03   3.670E-03   1.608E-03 /\n'
        '    CHEB /   8.599E-04  -1.758E-03  -7.502E-04   7.396E-07 /\n'
        'DUP\n\n\nEND\n\n')

    rxn_param_dct1 = get_rxn_param_dct(ckin_str1, 'cal/mole', 'moles')
    for params in rxn_param_dct1.values():
        cheb_dct = params.cheb
        tlim = cheb_dct['tlim']
        plim = cheb_dct['plim']
        alpha = cheb_dct['alpha']
        assert tlim == (500, 2000)
        assert plim == (0.03, 100)
        assert np.shape(alpha) == (6, 4)

    # Do the duplicates check
    rxn_param_dct2 = get_rxn_param_dct(ckin_str2, 'cal/mole', 'moles')
    for params in rxn_param_dct2.values():
        cheb_dct = params.cheb
        tlim = cheb_dct['tlim']
        plim = cheb_dct['plim']
        alpha = cheb_dct['alpha']
        cheb_dups = params.cheb_dups
        assert tlim == (500, 2000)
        assert plim == (0.03, 100)
        assert np.shape(alpha) == (6, 4)
        assert len(cheb_dups) == 1


def test_troe():
    """ Tests the Troe reader and writer
    """

    ckin_str1 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+12     1.500    50000\n'
        '    LOW  /              1.000E+12     1.500    50000  /\n'
        '    TROE /   1.500E+00   8.000E+03   1.000E+02   1.000E+03 /\n'
        '     AR/1.400/   N2/1.700/   \n\n\nEND\n\n')

    ckin_str2 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+12     1.500    50000\n'
        '    LOW  /              1.000E+12     1.500    50000  /\n'
        '    TROE /   1.500E+00   8.000E+03   1.000E+02   1.000E+03 /\n'
        '     AR/1.400/   N2/1.700/   \n'
        'DUP\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+12     1.500    50000\n'
        '    LOW  /              1.000E+12     1.500    50000  /\n'
        '    TROE /   1.500E+00   8.000E+03   1.000E+02   1.000E+03 /\n'
        '     AR/1.400/   N2/1.700/   \n'
        'DUP\n\n\nEND\n\n')

    rxn_param_dct1 = get_rxn_param_dct(ckin_str1, 'cal/mole', 'moles')
    for params in rxn_param_dct1.values():
        troe_dct = params.troe
        highp_arr = troe_dct['highp_arr']
        lowp_arr = troe_dct['lowp_arr']
        troe_params = troe_dct['troe_params']
        for arr_tuple in highp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])
        for arr_tuple in lowp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])
        assert np.allclose(troe_params, [1.5, 8e3, 1e2, 1e3])

    # Do the duplicates check
    rxn_param_dct2 = get_rxn_param_dct(ckin_str2, 'cal/mole', 'moles')
    for params in rxn_param_dct2.values():
        troe_dct = params.troe
        highp_arr = troe_dct['highp_arr']
        lowp_arr = troe_dct['lowp_arr']
        troe_params = troe_dct['troe_params']
        for arr_tuple in highp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])
        for arr_tuple in lowp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])
        assert np.allclose(troe_params, [1.5, 8e3, 1e2, 1e3])
        assert len(params.troe_dups) == 1


def test_lind():
    """ Tests the Lindemann reader and writer
    """

    ckin_str1 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+12     1.500    50000\n'
        '    LOW  /              1.000E+12     1.500    50000  /\n'
        '     AR/1.400/   N2/1.700/   \n\n\nEND\n\n')

    ckin_str2 = (
        'REACTIONS     CAL/MOLE     MOLES\n\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+12     1.500    50000\n'
        '    LOW  /              1.000E+12     1.500    50000  /\n'
        '     AR/1.400/   N2/1.700/   \n'
        'DUP\n'
        'H+O2(+N2)=OH+O(+N2)     1.000E+12     1.500    50000\n'
        '    LOW  /              1.000E+12     1.500    50000  /\n'
        '     AR/1.400/   N2/1.700/   \n'
        'DUP\n\n\nEND\n\n')

    rxn_param_dct1 = get_rxn_param_dct(ckin_str1, 'cal/mole', 'moles')
    for params in rxn_param_dct1.values():
        lind_dct = params.lind
        highp_arr = lind_dct['highp_arr']
        lowp_arr = lind_dct['lowp_arr']
        for arr_tuple in highp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])
        for arr_tuple in lowp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])

    # Do the duplicates check
    rxn_param_dct2 = get_rxn_param_dct(ckin_str2, 'cal/mole', 'moles')
    for params in rxn_param_dct2.values():
        lind_dct = params.lind
        highp_arr = lind_dct['highp_arr']
        lowp_arr = lind_dct['lowp_arr']
        for arr_tuple in highp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])
        for arr_tuple in lowp_arr:
            assert np.allclose(arr_tuple, [1e12, 1.5, 50000])
        assert len(params.lind_dups) == 1


def test_rxn_names():
    """ test mechanalyzer.parser.reaction.get_rxn_param_dct
        calls also
        get_rxn_strs
        get_rxn_name
        get_params
        fix_duplicates
    """
    NAMES_FROMDATA = {
        (('C2H3', 'O2'), ('CO', 'CH3O'), (None,)) ,
        (('C2H3', 'O2'), ('CO2', 'CH3'), (None,)) ,
        (('C2H3', 'O2'), ('CH2CHO', 'O'), (None,)) ,
        (('C2H3', 'O2'), ('C2H3OO',), (None,)) ,
        (('C2H3', 'O2'), ('C2H2', 'HO2'), (None,)) ,
        (('C2H3', 'O2'), ('CH2CO', 'OH'), (None,)) ,
        (('C2H3', 'O2'), ('CHOCHO', 'H'), (None,)) ,
        (('C2H3', 'O2'), ('CHCHO', 'OH'), (None,)) ,
        (('C2H3', 'O2'), ('CH2O', 'HCO'), (None,)) ,
    }
    ckin_str = ioformat.pathtools.read_file(DAT_PATH, 'rxn_block.dat')
    rxn_param_dct = get_rxn_param_dct(ckin_str, 'cal/mole', 'moles')
    assert all(name in rxn_param_dct.keys() for name in NAMES_FROMDATA)

def test_pes_dct():
    """ test mechanalyzer.parser.reaction.get_pes_dct with and w/o comments
        calls also
        get_rxn_strs
        get_rxn_name
        get_pes_info
    """
    RESULTS_PES_INFO = {('PES', 0, 0): ((0, (('C4H7ORvE4fmAB0',), ('C4H7O4H74fm1',), (None,))), 
                                        (1, (('C4H7ORvE4fmAB0',), ('C4H7O-kSV4fm',), (None,))), 
                                        (2, (('C4H7ORvE4fmAB0',), ('C4H6O-RvErx51', 'H-TcYTcY'), (None,))), 
                                        (3, (('C4H7ORvE4fmAA0',), ('C4H7O4H74fm0',), (None,))), 
                                        (4, (('C4H7ORvE4fmAA0',), ('C4H7O-kSV4fm',), (None,))), 
                                        (5, (('C4H7ORvE4fmAA0',), ('C4H6O-RvErx50', 'H-TcYTcY'), (None,))), 
                                        (6, (('C4H7O4H74fm0',), ('C3H4OALAD-Wv9FbZ', 'CH3'), (None,))), 
                                        (7, (('C4H7O4H74fm0',), ('C2H4OALD-UPQWKw', 'C2H3ALK-S58hH1'), (None,))), 
                                        (8, (('C4H7O-kSV4fm',), ('C2H4OALD-UPQWKw', 'C2H3ALK-S58hH1'), (None,)))), 
                        ('PES', 1, 0): ((0, (('C4H8ORvEsWvAA0', 'OH'), ('C4H7ORvE4fmAA0', 'H2O'), (None,))),), 
                        ('PES', 1, 1): ((1, (('C4H8ORvEsWvAB', 'OH'), ('C4H7ORvE4fmAB0', 'H2O'), (None,))),), 
                        ('PES', 2, 0): ((0, (('C4H8ORvEsWvAA0', 'HO2-S580KW'), ('C4H7ORvE4fmAA0', 'H2O2-S58pAY'), (None,))),), 
                        ('PES', 3, 0): ((0, (('C4H8ORvEsWvAA0', 'CH3'), ('C4H7ORvE4fmAA0', 'CH4'), (None,))),), 
                        ('PES', 4, 0): ((0, (('C4H8ORvEsWvAA0', 'CH3O2RO2-2LTcwB'), ('C4H7ORvE4fmAA0', 'CH4O2QOOH-2LTWKw'), (None,))),), 
                        ('PES', 5, 0): ((0, (('C4H8ORvEsWvAA0', 'CH3O-S58cwB'), ('C4H7ORvE4fmAA0', 'CH4O-S58WKw'), (None,))),), 
                        ('PES', 6, 0): ((0, (('C4H8ORvEsWvAA0', 'Cl'), ('C4H7ORvE4fmAA0', 'HCl'), (None,))),)}
    REACTION_BLOCK_PESUBPES_WDUP_INFO = {('PES', 0, 0): ((0, (('NC3H7',), ('IC3H7',), (None,))), 
    (1, (('C3H6', 'H'), ('NC3H7',), (None,))), 
    (2, (('C3H6', 'H'), ('IC3H7',), (None,))), 
    (3, (('C3H6', 'H'), ('C2H4', 'CH3'), (None,))), 
    (4, (('C2H4', 'CH3'), ('NC3H7',), (None,)))), 
    ('PES', 1, 0): ((0, (('C3H8', 'H'), ('IC3H7', 'H2'), (None,))),), 
    ('PES', 1, 1): ((1, (('C3H8', 'H'), ('NC3H7', 'H2'), (None,))),), 
    ('PES', 2, 0): ((0, (('C3H8', 'OH'), ('IC3H7', 'H2O'), (None,))),), 
    ('PES', 2, 1): ((1, (('C3H8', 'OH'), ('NC3H7', 'H2O'), (None,))),)}

    pes_info = get_pes_dct(REACTION_BLOCK_COMMENTS)
    assert pes_info == None
    pes_info = get_pes_dct(REACTION_BLOCK_OKCOMM)
    for key, val in pes_info.items():
        assert val == RESULTS_PES_INFO[key]
    pes_info = get_pes_dct(REACTION_BLOCK_PESUBPES_WDUP)
    for key, val in pes_info.items():
        assert val == REACTION_BLOCK_PESUBPES_WDUP_INFO[key]
if __name__ == '__main__':
    test_rxn_strs()
    test_cmts()
    test_rxn_osclass()
    test_nonint_rxn_products()
    test_arr()
    test_plog()
    test_plog2()
    test_cheb()
    test_troe()
    test_lind()
    test_rxn_names()
    test_pes_dct()
