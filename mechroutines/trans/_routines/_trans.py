"""
  CHEMKIN for ETRANS
"""

import automol
import autofile
import chemkin_io.writer
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
from mechroutines.models import etrans
from mechlib import filesys
from mechlib.amech_io import printer as ioprinter


def build_transport_file(tgt_queue,
                         spc_dct, thy_dct, pes_model_dct,
                         etrans_keyword_dct,
                         run_prefix, save_prefix):
    """ Write the chemkin string
    """

    # Get info from pes mod dct
    pes_mod = etrans_keyword_dct['kin_model']
    pes_model_dct_i = pes_model_dct[pes_mod]
    bath_name = pes_model_dct_i['energy_transfer']['bath']

    # Now obtain all the properties for the TGT+TGT interaction for CKIN
    # Uses simple combining rules for the LJ params
    trans_dct = {}
    for tgt_name in tgt_queue:

        print("\n+++++++++++++++++++++++++++++++++++++++++++++\n")
        print(f' Generating transport parameters for {tgt_name}')

        tgt_dct = spc_dct[tgt_name]
        bath_dct = spc_dct[bath_name]
        tgt_info = sinfo.from_dct(tgt_dct)

        # Build the filesystem objects needed for the TGT+BATH interaction
        tgt_cnf_inf, tb_etrans_inf, bb_etrans_inf = _etrans_fs(
            spc_dct, tgt_name, bath_name,
            thy_dct, etrans_keyword_dct,
            run_prefix, save_prefix)
        t_cnf_save_fs, t_min_cnf_locs = tgt_cnf_inf
        tb_etrans_fs, tb_etrans_locs = tb_etrans_inf
        bb_etrans_fs, bb_etrans_locs = bb_etrans_inf

        # Read the conformer filesystems
        if t_cnf_save_fs[-1].file.geometry.exists(t_min_cnf_locs):
            geo = t_cnf_save_fs[-1].file.geometry.read(t_min_cnf_locs)
            _path = t_cnf_save_fs[-1].file.geometry.path(t_min_cnf_locs)
            print(f'Reading geometry from path: {_path}')
        else:
            geo = None
        if t_cnf_save_fs[-1].file.dipole_moment.exists(t_min_cnf_locs):
            vec = t_cnf_save_fs[-1].file.dipole_moment.read(t_min_cnf_locs)
            dip_mom = automol.prop.total_dipole_moment(vec)
            _path = t_cnf_save_fs[-1].file.dipole_moment.path(t_min_cnf_locs)
            print(f'Reading dipole moment from path: {_path}')
        else:
            dip_mom = None
        if t_cnf_save_fs[-1].file.polarizability.exists(t_min_cnf_locs):
            tensor = t_cnf_save_fs[-1].file.polarizability.read(t_min_cnf_locs)
            polar = automol.prop.total_polarizability(tensor)
            _path = t_cnf_save_fs[-1].file.polarizability.path(t_min_cnf_locs)
            print(f'Reading polarizability from path: {_path}')
        else:
            polar = None

        # Obtain LJ sigma, epsilon from filesystem or estimation
        print('Determining Target+Bath and Bath+Bath LJ parameters...')
        tb_sig, tb_eps, bb_sig, bb_eps = _lj_params(
            tb_etrans_fs, tb_etrans_locs,
            bb_etrans_fs, bb_etrans_locs,
            tgt_dct, bath_dct,
            pes_model_dct_i)

        # Use the combining rules for sigma and epsilon
        print('Determining Target+Target LJ parameters '
              'from combining rules...')
        tt_eps = automol.etrans.combine.epsilon(tb_eps, bb_eps)
        tt_sig = automol.etrans.combine.sigma(tb_sig, bb_sig)

        # Get the Z_ROT number
        print('Determining Z ROT at 298 K...')
        zrot = automol.etrans.estimate.rotational_relaxation_number(
            tgt_info[0])

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
    ioprinter.info_message(
        'Collating data into CHEMKIN transport file string:', newline=1)
    transport_str = chemkin_io.writer.transport.properties(trans_dct)
    ioprinter.debug_message(transport_str)
    ioprinter.obj('vspace')
    ioprinter.obj('line_dash')

    return transport_str


def _lj_params(tb_etrans_fs, tb_etrans_locs,
               bb_etrans_fs, bb_etrans_locs,
               tgt_dct, bath_dct,
               pes_mod_dct_i):
    """ Obtain the energy transfer parameters by reading the filesystem
        or estimating it with the Jasper formulae.
    """

    if (
        tb_etrans_fs[-1].file.lennard_jones_epsilon.exists(tb_etrans_locs) and
        bb_etrans_fs[-1].file.lennard_jones_epsilon.exists(bb_etrans_locs) and
        tb_etrans_fs[-1].file.lennard_jones_sigma.exists(tb_etrans_locs) and
        bb_etrans_fs[-1].file.lennard_jones_sigma.exists(bb_etrans_locs)
    ):
        print('Reading lennard_jones_sigma, lennard_jones_epsilon from filesystem')
        tb_eps = tb_etrans_fs[-1].file.lennard_jones_epsilon.read(tb_etrans_locs)
        bb_eps = bb_etrans_fs[-1].file.lennard_jones_epsilon.read(bb_etrans_locs)
        tb_sig = tb_etrans_fs[-1].file.lennard_jones_sigma.read(tb_etrans_locs)
        bb_sig = bb_etrans_fs[-1].file.lennard_jones_sigma.read(bb_etrans_locs)
    else:
        t_etrans_dct = etrans.etrans_dct_for_species(tgt_dct, pes_mod_dct_i)
        b_etrans_dct = etrans.etrans_dct_for_species(bath_dct, pes_mod_dct_i)
        tgt_info = sinfo.from_dct(tgt_dct)
        bath_info = sinfo.from_dct(bath_dct)
        tb_sig, tb_eps = etrans.lj_params(tgt_info, bath_info, t_etrans_dct)
        bb_sig, bb_eps = etrans.lj_params(bath_info, bath_info, b_etrans_dct)

    return bb_sig, bb_eps, tb_sig, tb_eps


def _etrans_fs(spc_dct, tgt_name, bath_name,
               thy_dct, etrans_keyword_dct,
               run_prefix, save_prefix):
    """ Build the energy transfer filesys

        Need to get the cnf_fs for the target
        and the etrans_fs for the target-bath, bath-bath
    """

    method_dct = thy_dct.get(etrans_keyword_dct['inplvl'])
    thy_info = tinfo.from_dct(method_dct)
    bath_info = sinfo.from_dct(spc_dct[bath_name])

    cnf_fs_lst, cnf_locs_lst = [], []
    etrans_fs_lst, etrans_locs_lst = [], []
    for name in (tgt_name, bath_name):
        spc_dct_i = spc_dct[name]
        spc_info = sinfo.from_dct(spc_dct_i)
        mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)

        # Build the conformer filesystem objects
        cnf_sort_info_lst = _sort_info_lst(
            etrans_keyword_dct['sort'], thy_dct, spc_info)

        _, cnf_save_fs = filesys.build_fs(
            run_prefix, save_prefix, 'CONFORMER',
            spc_locs=spc_info,
            thy_locs=mod_thy_info[1:])

        cnf_rng_info = filesys.mincnf.conformer_locators(
            cnf_save_fs, mod_thy_info,
            cnf_range=etrans_keyword_dct['cnf_range'],
            sort_info_lst=cnf_sort_info_lst,
            hbond_cutoffs=spc_dct_i['hbond_cutoffs'],
            print_enes=True,
            nprocs=1)
        min_cnf_locs = cnf_rng_info[0][0]

        # Build the energy transfer filesystem objects
        cnf_path = cnf_save_fs[-1].path(min_cnf_locs)
        etrans_fs = autofile.fs.energy_transfer(cnf_path)

        # Build the energy transfer locs
        lj_info = sinfo.combine(spc_info, bath_info)
        mod_lj_thy_info = tinfo.modify_orb_label(thy_info, lj_info)

        etrans_locs = bath_info + mod_lj_thy_info[1:]

        # Add to the lists
        cnf_fs_lst.append(cnf_save_fs)
        cnf_locs_lst.append(min_cnf_locs)
        etrans_fs_lst.append(etrans_fs)
        etrans_locs_lst.append(etrans_locs)

    return (
        (cnf_fs_lst[0], cnf_locs_lst[0]),                # target
        (etrans_fs_lst[0], etrans_locs_lst[0]),          # target-bath
        (etrans_fs_lst[1], etrans_locs_lst[1]),          # bath-bath
    )


def _sort_info_lst(sort_str, thy_dct, spc_info):
    """ Return the levels to sort conformers by if zpve or sp
        levels were assigned in input

        if we ask for zpe(lvl_wbs),sp(lvl_b2t),gibbs(700)
        out sort_info_lst will be [('gaussian', 'wb97xd', '6-31*', 'RU'),
        ('gaussian', 'b2plypd3', 'cc-pvtz', 'RU'), None, None, 700.]
    """
    sort_lvls = [None, None, None, None, None]
    sort_typ_lst = ['freqs', 'sp', 'enthalpy', 'entropy', 'gibbs']
    if sort_str is not None:
        for sort_param in sort_str.split(','):
            idx = None
            for typ_idx, typ_str in enumerate(sort_typ_lst):
                if typ_str in sort_param:
                    lvl_key = sort_str.split(typ_str + '(')[1].split(')')[0]
                    idx = typ_idx
            if idx is not None:
                if idx < 2:
                    method_dct = thy_dct.get(lvl_key)
                    if method_dct is None:
                        ioprinter.warning_message(
                            f'no {lvl_key} in theory.dat, '
                            f'not using {sort_typ_lst[idx]} in sorting')
                        continue
                    thy_info = tinfo.from_dct(method_dct)
                    mod_thy_info = tinfo.modify_orb_label(thy_info, spc_info)
                    sort_lvls[idx] = mod_thy_info
                else:
                    sort_lvls[idx] = float(lvl_key)
    return sort_lvls
