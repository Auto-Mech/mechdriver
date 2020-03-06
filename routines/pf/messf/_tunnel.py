"""
Interface to MESS and projrot to set-up tunneling blocks
"""

<<<<<<< HEAD
import os
import mess_io
import projrot_io
import autofile
from lib.runner import script
from lib.filesystem import orb as fsorb
from routines.es.variational import irc

# rxn path coords from g09 outpt
# rxn path enes from g09 outpt


def write_mess_eckart_str(ts_ene, reac_ene, prod_ene, imag_freq):
    """ Write the Eckart tunneling string for MESS'
    """

    # Set the depth of the wells from the transition state
    ts_reac_barr = ts_ene - reac_ene
    ts_prod_barr = ts_ene - prod_ene
    if ts_reac_barr < 0.:
        ts_reac_barr = 0.1
    if ts_prod_barr < 0.:
        ts_prod_barr = 0.1

    # Write the MESS string
    tunnel_str = mess_io.writer.tunnel_eckart(
        imag_freq, ts_reac_barr, ts_prod_barr)

    return tunnel_str


def write_mess_sct_str(ts_dct, pf_levels, save_path,
                       imag_freq, tunnel_file,
                       cutoff_energy=2500, coord_proj='cartesian'):
    """ Write the Eckart tunneling string for MESS
    """
    tunnel_str = mess_io.writer.tunnel_sct(
        imag_freq, tunnel_file, cutoff_energy=cutoff_energy)
    trans_coeff_str = build_trans_coeff_file(
        ts_dct, pf_levels, save_path, coord_proj=coord_proj)

    return tunnel_str, trans_coeff_str


def build_trans_coeff_file(ts_dct, pf_levels,
                           save_path, coord_proj='cartesian'):
    """ Collate the IRC data and build
        the transmission coefficient file with ProjRot
    """

    # Loop over the IRC filesys and read the info
    [geo_thy_level, ene_thy_level, _, _, _, _] = pf_levels

    print('ts_dct')
    print(ts_dct)
    irc_idxs = ts_dct['irc_idxs']
    ts_info = ['', ts_dct['chg'], ts_dct['mul']]

    # Build the TS scan file system
    scn_save_fs, _, _, _ = irc.ts_scn_fs(
        ts_dct, ts_info, geo_thy_level)

    # Set the distance name for the reaction coordinate
    dist_name = 'RC'

    # Loop over all the points of the irc and build MESS strings
    rpath_coords = []
    geoms, grads, hessians, enes = [], [], [], []
    for idx in irc_idxs:

        # Set the filesystem locators for each grid point
        locs = [[dist_name], [idx]]
        print(scn_save_fs[-1].path(locs))

        # Get geometries, gradients, hessians, enes
        if scn_save_fs[-1].file.geometry.exists(locs):
            geoms.append(scn_save_fs[-1].file.geometry.read(locs))
        else:
            print('no geom')
            continue
        if scn_save_fs[-1].file.energy.exists(locs):
            if ene_thy_level == geo_thy_level:
                enes.append(scn_save_fs[-1].file.energy.read(locs))
                print('geo ene read')
            else:
                scn_save_path = scn_save_fs[-1].path(locs)
                sp_save_fs = autofile.fs.single_point(scn_save_path)
                orb_restr = fsorb.orbital_restriction(ts_info, ene_thy_level)
                sp_level = ene_thy_level[0:3]
                sp_level.append(orb_restr)
                if sp_save_fs[-1].file.energy.exists(sp_level[1:4]):
                    enes.append(sp_save_fs[-1].file.energy.read(sp_level[1:4]))
                else:
                    print('no energy')
                    continue
        else:
            print('no energy')
            continue
        if scn_save_fs[-1].file.gradient.exists(locs):
            grads.append(scn_save_fs[-1].file.gradient.read(locs))
        else:
            print('no gradient')
            continue
        if scn_save_fs[-1].file.hessian.exists(locs):
            hessians.append(scn_save_fs[-1].file.hessian.read(locs))
        else:
            print('no hessian')
            continue

        # Set the saddle index and energy
        if idx == 0.00:
            saddle_idx = irc_idxs.index(idx) + 1
            saddle_ene = enes[-1]

        # append coords for now
        rpath_coords.append(0.01)

    # Scale the energies
    rpath_enes = [ene - saddle_ene for ene in enes]

    # Write the input and coord-energy string for the ProjRot input
    projrot_inp_str = projrot_io.writer.rpht_input(
        geoms, grads, hessians,
        saddle_idx=saddle_idx,
        rotors_str='', coord_proj=coord_proj)
    projrot_en_str = projrot_io.writer.rpht_path_coord_en(
        rpath_coords, rpath_enes,
        bnd1=(), bnd2=())

    print('SCT test')
    print(projrot_inp_str)
    print(projrot_en_str)

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
    with open(proj_file_path, 'w') as proj_file:
        proj_file.write(projrot_en_str)
    script.run_script(script.PROJROT, path)

    # Read the transmission coefficient file
    with open(os.path.join(path, 'imactint.txt')) as tc_file:
        trans_coeff_str = tc_file.read()

    return trans_coeff_str
=======
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




>>>>>>> 9bf6c45022eda7c7b6875f3290b6564dd75cf7b9
