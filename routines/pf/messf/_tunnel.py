"""
Interface to MESS and projrot to set-up tunneling blocks
"""

import os
import mess_io
import projrot_io
import autofile
from lib import filesys
from runners import run_script
from runners import DEFAULT_SCRIPT_DCT


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
                orb_restr = filesys.inf.orbital_restriction(ts_info, ene_thy_level)
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
    run_script(DEFAULT_SCRIPT_DCT['projrot'], path)

    # Read the transmission coefficient file
    with open(os.path.join(path, 'imactint.txt')) as tc_file:
        trans_coeff_str = tc_file.read()

    return trans_coeff_str


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
