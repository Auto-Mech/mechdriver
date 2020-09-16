"""
  CHEMKIN for ETRANS
"""

import automol
import chemkin_io
from lib import filesys


def collate_properties(spc_queue, bath_name,
                       spc_dct, thy_dct, etrans_keyword_dct,
                       save_prefix):
    """ Write the chemkin string
    """

    # Get the base theory info obj
    thy_info = filesys.inf.get_es_info(
        etrans_keyword_dct['runlvl'], thy_dct)

    # Build the bath info object
    bath_info = filesys.inf.get_spc_info(spc_dct[bath_name])

    names, geoms = [], []
    dipole_moments, polarizabilities = [], []
    epsilons, sigmas = [], []
    for spc_name, _ in spc_queue:

        # Get the info for the target combined spc info objects
        tgt_dct = spc_dct[spc_name]
        tgt_info = filesys.inf.get_spc_info(tgt_dct)
        lj_info = filesys.inf.combine_spc_info(tgt_info, bath_info)

        # Build the modified thy objs
        mod_tgt_thy_info = filesys.inf.modify_orb_restrict(tgt_info, thy_info)
        mod_lj_thy_info = filesys.inf.modify_orb_restrict(lj_info, thy_info)

        # Build the conformer filesystem objects
        _, thy_save_path = filesys.build.spc_thy_fs_from_root(
            save_prefix, tgt_info, mod_tgt_thy_info)
        cnf_save_fs, cnf_save_locs = filesys.build.cnf_fs_from_prefix(
            thy_save_path, mod_tgt_thy_info, cnf='min')
        cnf_save_paths = filesys.build.cnf_paths_from_locs(
            cnf_save_fs, cnf_save_locs)

        # Build the energy transfer filesystem objects
        etrans_fs, etrans_locs = filesys.build.etrans_fs_from_prefix(
            cnf_save_paths[0], bath_info, mod_lj_thy_info)

        # Read the conformer filesystems
        if cnf_save_fs[-1].file.geometry.exists(cnf_save_locs):
            geom = cnf_save_fs[-1].file.geometry.read(cnf_save_locs)
        else:
            geom = None
        if cnf_save_fs[-1].file.dipole_moment.exists(cnf_save_locs):
            vec = cnf_save_fs[-1].file.dipole_moment.read(cnf_save_locs)
            dip_mom = automol.prop.total_dipole_moment(vec)
        else:
            dip_mom = None
        if cnf_save_fs[-1].file.polarizability.exists(cnf_save_locs):
            tensor = cnf_save_fs[-1].file.polarizability.read(cnf_save_locs)
            polar = automol.prop.total_polarizability(tensor)
        else:
            polar = None

        # Read the energy transfer filesystems
        if etrans_fs[-1].file.epsilon.exists(etrans_locs):
            eps = etrans_fs[-1].file.epsilon.read(etrans_locs)
        else:
            eps = None
        if etrans_fs[-1].file.sigma.exists(etrans_locs):
            sig = etrans_fs[-1].file.sigma.read(etrans_locs)
        else:
            sig = None

        # Append info to list
        names.append(spc_name)
        geoms.append(geom)
        dipole_moments.append(dip_mom)
        polarizabilities.append(polar)
        epsilons.append(eps)
        sigmas.append(sig)

    # Write the string with all of the transport properties
    transport_str = chemkin_io.writer.transport.properties(
        names, geoms, epsilons, sigmas,
        dipole_moments, polarizabilities, z_rots=None)
    print('\ntransport_str')
    print(transport_str)

    return transport_str
