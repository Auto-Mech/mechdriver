"""
  CHEMKIN for ETRANS
"""

import automol
import autofile
import chemkin_io.writer
from mechroutines.models import etrans
from mechlib import filesys
from mechlib.amech_io import printer as ioprinter


def build_transport_file(tgt_queue,
                         spc_dct, thy_dct, etrans_keyword_dct,
                         save_prefix):
    """ Write the chemkin string
    """

    bath_name = etrans_keyword_dct['bath']

    # Read the epsilon and sigma params for the BATH+BATH interaction
    _, _, bb_etrans_fs, bb_etrans_locs = _etrans_fs(
        spc_dct, bath_name, bath_name,
        thy_dct, etrans_keyword_dct,
        save_prefix)

    # Now obtain all the properties for the TGT+TGT interaction for CKIN
    # Uses simple combining rules for the LJ params
    trans_dct = {}
    for tgt_name, _ in tgt_queue:

        # Build the filesystem objects needed for the TGT+BATH interaction
        cnf_save_fs, min_cnf_locs, tt_etrans_fs, tt_etrans_locs = _etrans_fs(
            spc_dct, tgt_name, bath_name,
            thy_dct, etrans_keyword_dct,
            save_prefix)

        # Read the conformer filesystems
        if cnf_save_fs[-1].file.geometry.exists(min_cnf_locs):
            geo = cnf_save_fs[-1].file.geometry.read(min_cnf_locs)
        else:
            geo = None
        if cnf_save_fs[-1].file.dipole_moment.exists(min_cnf_locs):
            vec = cnf_save_fs[-1].file.dipole_moment.read(min_cnf_locs)
            dip_mom = automol.prop.total_dipole_moment(vec)
        else:
            dip_mom = None
        if cnf_save_fs[-1].file.polarizability.exists(min_cnf_locs):
            tensor = cnf_save_fs[-1].file.polarizability.read(min_cnf_locs)
            polar = automol.prop.total_polarizability(tensor)
        else:
            polar = None

        # Obtain LJ sigma, epsilon from filesystem or estimation
        bb_sig, bb_eps, tb_sig, tb_eps = _lj_params(
            bb_etrans_fs, tt_etrans_fs, bb_etrans_locs, tt_etrans_locs,
            tgt_info, bath_info)

        # Use the combining rules for sigma and epsilon
        tt_eps = automol.etrans.combine.epsilon(tb_eps, bb_eps)
        tt_sig = automol.etrans.combine.sigma(tb_sig, bb_sig)

        # Get the Z_ROT number
        zrot = automol.etrans.rotational_relaxation_number(tgt_info)

        # Build dct and append it ot overall dictionary Append info to list
        dct = {
            'geo': geo,
            'dipole_moment': dip_mom,
            'polarizability': polar,
            'epsilon': tt_eps,
            'sigma': tt_sig,
            'zrot': zrot
        }
        trans_dct.update({tgt_name: dct})

    # Write the string with all of the transport properties
    transport_str = chemkin_io.writer.transport.properties(trans_dct)
    ioprinter.debug_message('transport_str\n', transport_str, newline=1)
    ioprinter.obj('vspace')
    ioprinter.obj('line_dash')
    ioprinter.info_message('Writing the CHEMKIN transport file', newline=1)

    return transport_str


def _lj_params(bb_etrans_fs, tb_etrans_fs, bb_etrans_locs, tb_etrans_locs,
               tgt_info, bath_info,
               pes_model_dct_i):
    """ Obtain the energy transfer parameters by reading the filesystem
        or estimating it with the Jasper formulae.
    """

    etransfer = pes_model_dct_i['glob_etransfer']

    if (
        bb_etrans_fs[-1].file.epsilon.exists(bb_etrans_locs) and
        tb_etrans_fs[-1].file.epsilon.exists(tb_etrans_locs) and
        bb_etrans_fs[-1].file.sigma.exists(bb_etrans_locs) and
        tb_etrans_fs[-1].file.sigma.exists(tb_etrans_locs)
    ):
        bb_eps = bb_etrans_fs[-1].file.epsilon.read(bb_etrans_locs)
        tb_eps = tb_etrans_fs[-1].file.epsilon.read(tb_etrans_locs)
        bb_sig = bb_etrans_fs[-1].file.sigma.read(bb_etrans_locs)
        tb_sig = tb_etrans_fs[-1].file.sigma.read(tb_etrans_locs)
    else:
        etrans_dct = pes_model_dct_i['glob_etransfer']
        tb_sig, tb_eps, _, _ = etrans.lj_params(
            tgt_info, bath_info, etrans_dct)
        bb_sig, bb_eps, _, _ = etrans.lj_params(
            tgt_info, bath_info, etrans_dct)

    return bb_sig, bb_eps, tb_sig, tb_eps


def _etrans_fs(spc_dct, tgt_name, bath_name,
               thy_dct, etrans_keyword_dct,
               save_prefix):
    """ Build the energy transfer filesys
    """

    # Get the base theory info obj
    thy_info = filesys.inf.get_es_info(
        etrans_keyword_dct['runlvl'], thy_dct)

    min_locs = ()
    spc_infos = ()
    for name in (tgt_name, bath_name):
        spc_dct_i = spc_dct[name]
        spc_info = filesys.inf.get_spc_info(spc_dct_i)
        mod_thy_info = filesys.inf.modify_orb_restrict(spc_info, thy_info)

        # Build the conformer filesystem objects
        cnf_range = etrans_keyword_dct['cnf_range']
        hbond_cutoffs = spc_dct_i['hbond_cutoffs']
        cnf_sort_info_lst = _sort_info_lst(
            etrans_keyword_dct['sort'], thy_dct, spc_info)

        _, cnf_save_fs = filesys.build_fs(
            run_prefix, save_prefix, 'CONFORMER',
            spc_locs=spc_info, thy_locs=mod_thy_info[1:])

        cnf_rng_info = filesys.mincnf.conformer_locators(
            cnf_save_fs, mod_thy_info,
            cnf_range=etrans_keyword_dct['cnf_range'],
            sort_info_lst=cnf_sort_info_lst,
            hbond_cutoffs=spc_dct_i['hbond_cutoffs'],
            print_enes=True,
            nprocs=1)
        min_cnf_locs, min_cnf_path = cnf_rng_info

    # Build the energy transfer filesystem objects
    lj_info = filesys.inf.combine_spc_info(*spc_infos)
    mod_lj_thy_info = filesys.inf.modify_orb_restrict(lj_info, thy_info)
    _, etrans_save_fs = filesys.build_fs(
        run_prefix, save_prefix, 'ENERGY TRANSFER',
        spc_locs=spc_info,
        thy_locs=mod_lj_thy_info[1:],
        cnf_locs=min_cnf_locs)

    return cnf_save_fs, min_cnf_locs, etrans_fs, etrans_locs
