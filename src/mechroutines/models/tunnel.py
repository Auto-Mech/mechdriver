"""
NEW: Interface to MESS and projrot to set-up tunneling blocks
"""

import mess_io
from mechroutines.models.typ import treat_tunnel


def write_mess_tunnel_str(ts_inf_dct, chnl_enes,
                          ts_model, ts_class, ts_idx, unstable_chnl=False):
    """ Write a `Tunneling` section string for a transition state for a MESS
        input file for the requested tunneling model.

        Transition state data from the SAVE filesystem is processed to produce
        the values required to produce the MESS string for the appropriate
        model. When needed, external codes such as ProjRot or PolyRate may
        be called to produce this required data for the MESS strings.

        :param ts_inf_dct: collated transition state info from SAVE filesys
        :type ts_inf_dct: dict[str:__]
        :param chnl_enes: energies for channel, relative to PES reference
        :type chnl_enes: dict[str:float]
        :param ts_model: transition state treatment model from spc model dct
        :type ts_model: dict[str: str]
        :param ts_idx: idx along eaction path denoting position of saddle point
        :type ts_idx: int
        :param unstable_chnl: denotes if unstable species in reacs/prods
        :type unstable_chnl: bool
    """

    tunnel_str, sct_dat = '', {}
    if treat_tunnel(ts_model, ts_class, ts_inf_dct):
        tunnel_model = ts_model['tunnel']
        if tunnel_model == 'eckart':
            tunnel_str = write_mess_eckart_str(
                chnl_enes, ts_inf_dct.get('imag', None),
                ts_idx=ts_idx, symm_barrier=unstable_chnl)
        # elif tunnel_model == 'sct':
        #     sct_dat_name = tsname + '_sct.dat'
        #     path = 'cat'
        #     tunnel_str, sct_str = write_mess_sct_str(
        #         spc_dct[tsname], pf_levels, path,
        #         imag, sct_dat_name,
        #         cutoff_energy=2500, coord_proj='cartesian')

        #     sct_dat = {sct_dat_name: sct_str}

    return tunnel_str, sct_dat


def write_mess_eckart_str(chnl_enes, imag_freq, ts_idx=0, symm_barrier=False):
    """ Write a Eckart model `Tunneling` section string for a transition state
        for a MESS input file.

        Function simply processes the energies of the reactants and products
        to calculate the well-depths relative to the saddle point of the
        reaction coordinate. When needed, the depths are set to be similar
        for unstable reactants or products. For cases involving very small
        well-depths are calculated, the depths are set to be 0.1 to avoid
        unphysical master equation simulations..

        :param chnl_enes: energies for channel, relative to PES reference
        :type chnl_enes: dict[str:float]
        :param imag_freq: imaginary frequency of reaction coordinate [cm-1]
        :type imag_freq: float
        :param ts_idx: idx along eaction path denoting position of saddle point
        :type ts_idx: int
        :param symm_barrier: set the well-depths of Eckart model to be similar
        :type symm_barrier: bool
        :rtype: str
    """

    # Get the energies from the enes dct
    ts_ene = chnl_enes['ts'][ts_idx]

    if chnl_enes.get('fake_vdwr', None) is not None:
        reac_ene = chnl_enes['fake_vdwr']
    else:
        reac_ene = chnl_enes['reacs']

    if symm_barrier:
        prod_ene = reac_ene
    else:
        if chnl_enes.get('fake_vdwp', None) is not None:
            prod_ene = chnl_enes['fake_vdwp']
        else:
            prod_ene = chnl_enes['prods']

    # Set the depth of the wells from the transition state
    ts_reac_barr = ts_ene - reac_ene
    ts_prod_barr = ts_ene - prod_ene
    ts_reac_barr = ts_reac_barr if ts_reac_barr > 0.1 else 0.1
    ts_prod_barr = ts_prod_barr if ts_prod_barr > 0.1 else 0.1

    # Write the MESS string
    tunnel_str = mess_io.writer.tunnel_eckart(
        imag_freq, ts_reac_barr, ts_prod_barr)

    return tunnel_str


# def write_mess_sct_strs(ts_inf_dct, save_path,
#                         imag_freq, tunnel_file, sct_prog='projrot'):
#     """ Write the Eckart tunneling string for MESS
#     """
#     # Write the tunneling section of TS for MESS input file
#     tunnel_str = mess_io.writer.tunnel_sct(
#         imag_freq, tunnel_file, cutoff_energy=2500)
#     # Write the auxiliary file with the tunneling action
#     if sct_prog == 'projrot':
#         tc_str = _projrot_tunn_action_str(ts_inf_dct, save_path)
#     elif sct_prog == 'polyrate':
#         tc_str = _polyrate_tunn_action_str(ts_inf_dct, save_path)
#     return tunnel_str, tc_str
# def _projrot_tunn_action_str(ts_inf_dct, run_prefix):
#     """ Collate the IRC data and build
#         the transmission coefficient file with ProjRot
#     """
#     # Loop over all the points of the irc and build MESS strings
#     sadpt_idx = ts_inf_dct['rpath']['ts_idx'] + 1
#     sadpt_ene = ts_inf_dct['rpath'][ts_idx]['ene_chnlvl']
#     # Obtain the energies and coords
#     rpath_enes, rpath_coords = [], []
#     for _, dct in enumerate(inf_dct['rpath']):
#         rpath_enes.append(dct['ene_chnlvl'] - sadpt_ene)
#         rpath_coords.append(dct['rvals'])
#     fml_str = automol.geom.formula_string(harm_geo)
#     sct_path = job_path(run_prefix, 'PROJROT', 'SCT',
#                         fml_str, print_path=True)
#     tc_str = autorun.projrot.small_curvature_tunneling(
#         script_str, run_dir, geoms, grads, hessians,
#         rpath_coords, rpath_enes, sadpt_idx,
#         rotors_str='')
#     return tc_str
# def _polyrate_tunn_action_str():
#     """ calculate s(E) using polyrate
#     """
#     raise NotImplementedError
