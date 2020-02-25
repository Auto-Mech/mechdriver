"""
Interface to MESS and projrot to set-up tunneling blocks
"""

RXN_PATH_COORDS = [
    0.30596, 0.27536, 0.24477, 0.21417, 0.18358,
    0.15298, 0.12238, 0.09179, 0.06119, 0.03060,
    0.00000,
    -0.03060, -0.06119, -0.09178, -0.12237, -0.15296,
    -0.18354, -0.21412, -0.24471, -0.27530, -0.30589]
RXN_PATH_ENERGIES = [
    -0.0046300, 0.0038600, 0.0031100, 0.0024100, 0.0017800,
    -0.0012300, 0.0007800, 0.0004300, 0.0001900, 0.0000400,
    0.0000000,
    -0.0000400, -0.0001500, -0.0003300, -0.0005400, -0.0007800,
    -0.0010300, -0.0013000, -0.0015600, -0.0018300, -0.0020900]
RCT_DISTS = [
    1.31974, 1.30370, 1.28780, 1.27204, 1.25643,
    1.24095, 1.22561, 1.21045, 1.19547, 1.18071,
    1.16640,
    1.15214, 1.13850, 1.12545, 1.11313, 1.10176,
    1.09156, 1.08265, 1.07505, 1.06868, 1.06336]
PRD_DISTS = [
    1.08945, 1.10474, 1.12025, 1.13593, 1.15176,
    1.16770, 1.18375, 1.19987, 1.21609, 1.23238,
    1.24880,
    1.26522, 1.28171, 1.29821, 1.31466, 1.33097,
    1.34701, 1.36265, 1.37782, 1.39250, 1.40671]


def write_rpht_sct_input(coord_proj='cartesian'):
    """ Read the struct info for each point alon
    """

    # Write the string for the ProjRot input
    inp_str = projrot_io.writer.rpht_input(
        GEOMS, GRADS, HESSES,
        saddle_idx=11, rotors_str='', coord_proj=CART_PROJ)
      
    # Print the string

    return inp_str


def write_rpht_sct_coord_en():
    """ test projrot_io.writer.rpht_path_coord_en
    """
    # Write the string withoutp bnd1 and bnd2 vals
    en_str = projrot_io.writer.rpht_path_coord_en(
        RXN_PATH_COORDS, RXN_PATH_ENERGIES,
        bnd1=None, bnd2=None)


    return en_str

    # Print the string
    print(en_str)

    # Write the string with bnd1 and bnd2 vals
    en_str = projrot_io.writer.rpht_path_coord_en(
        RXN_PATH_COORDS, RXN_PATH_ENERGIES,
        bnd1=RCT_DISTS, bnd2=PRD_DISTS)

    # Print the string
    print(en_str)




