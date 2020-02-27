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


def build_trans_coeff_file():
    """ Collate the IRC data and build 
        the transmission coefficient file with ProjRot
    """

    # Loop over the IRC filesys and read the info


    # Write the input and coord-energy string for the ProjRot input
    inp_str = projrot_io.writer.rpht_input(
        geoms, grads, hesses,
        saddle_idx=11, rotors_str='', coord_proj=cart_proj)
    en_str = projrot_io.writer.rpht_path_coord_en(
        rxn_path_coords, rxn_path_energies,
        bnd1=none, bnd2=none)

    # Write the ProjRot input files and run the code 
    bld_locs = ['PROJROT', 0]
    bld_save_fs = autofile.fs.build(save_path)
    bld_save_fs[-1].create(bld_locs)
    path = bld_save_fs[-1].path(bld_locs)
    print('Build Path for ProjRot calls')
    print(path)
    proj_file_path = os.path.join(path, 'RPHt_input_data.dat')
    with open(proj_file_path, 'w') as proj_file:
        proj_file.write(projrot_inp_str)
    script.run_script(script.PROJROT, path)

    return en_str

