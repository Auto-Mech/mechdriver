"""
NEW: Interface to MESS and projrot to set-up tunneling blocks
"""

import mess_io


def write_mess_eckart_str(chnl_enes, imag_freq, ts_idx=0):
    """ Write the Eckart tunneling string for MESS'
    """

    # Get the energies from the enes dct
    ts_ene = chnl_enes['ts'][ts_idx]
    # print('chnl_enes test:', chnl_enes)

    if chnl_enes.get('fake_vdwr', None) is not None:
        reac_ene = chnl_enes['fake_vdwr']
    else:
        reac_ene = chnl_enes['reacs']

    if chnl_enes.get('fake_vdwp', None) is not None:
        prod_ene = chnl_enes['fake_vdwp']
    else:
        prod_ene = chnl_enes['prods']

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


def write_mess_sct_strs(ts_dct, pf_levels, save_path,
                        imag_freq, tunnel_file,
                        cutoff_energy=2500, coord_proj='cartesian'):
    """ Write the Eckart tunneling string for MESS
    """
    tunnel_str = mess_io.writer.tunnel_sct(
        imag_freq, tunnel_file, cutoff_energy=cutoff_energy)
    trans_coeff_str = build_trans_coeff_file(
        ts_dct, pf_levels, save_path, coord_proj=coord_proj)

    return tunnel_str, trans_coeff_str


def build_trans_coeff_file(ts_inf_dct, run_path, coord_proj='cartesian'):
    """ Collate the IRC data and build
        the transmission coefficient file with ProjRot
    """

    # Loop over all the points of the irc and build MESS strings
    sadpt_idx = ts_inf_dct['rpath']['ts_idx'] + 1
    sadpt_ene = ts_inf_dct['rpath'][ts_idx]['ene_chnlvl']

    # Obtain the energies and coords
    rpath_enes, rpath_coords = [], []
    for idx, dct in enumerate(inf_dct['rpath']):
        rpath_enes.append(dct['ene_chnlvl'] - sadpt_ene)
        rpath_coords.append(dct['rvals'])

    # Write the input and coord-energy string for the ProjRot input
    projrot_inp_str = projrot_io.writer.rpht_input(
        geoms, grads, hessians,
        saddle_idx=sadpt_idx,
        rotors_str='',
        coord_proj=coord_proj)
    projrot_en_str = projrot_io.writer.rpht_path_coord_en(
        rpath_coords, rpath_enes,
        bnd1=(), bnd2=())

    print('SCT test')
    print(projrot_inp_str)
    print(projrot_en_str)

    # Write the ProjRot input files and run the code
    bld_locs = ['PROJROT', 0]
    bld_save_fs = autofile.fs.build(run_path)
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
