"""
  CHEMKIN for ETRANS
"""

import automol
import chemkin_io
from lib import filesys


# Put in tsk file or new lib

def collate_properties(spc_queue, spc_name, bath_name,
                       spc_dct, pf_levels,
                       thy_dct, etrans_keyword_dct,
                       run_prefix, save_prefix):
    """ Write the chemkin string
    """

    thy_info = filesys.inf.get_es_info(
        etrans_keyword_dct['inplvl'], thy_dct)
    # Need to get a combined spin for the modification
    mod_thy_info = filesys.inf.modify_orb_restrict(spc_info, thy_info)

    names, geoms = [], []
    dipole_moments, polarizabilities = [], []
    epsilons, sigmas = [], []
    for spc_name, _ in spc_queue:

        tgt_dct = spc_dct[spc_name]
        bath_dct = spc_dct[bath_name]

        bath_info = get_spc_info(spc_dct[bath_name])



        # Build the filesystem objects
        pf_filesystems = pf_filesys(
            tgt_dct, pf_levels,
            run_prefix, save_prefix, saddle=False)
        [cnf_fs, _, cnf_locs, _, _] = pf_filesystems['harm']

        etrans_fs, etrans_locs = _etrans_fs(
            cnf_fs, cnf_locs, bath_dct, thy_info)

        # Read the conformer filesystems
        if cnf_fs[-1].file.geometry.exist(cnf_locs):
            geom = cnf_fs[-1].file.geometry.read(cnf_locs)
        else:
            geom = None
        if cnf_fs[-1].file.dipole_moment.exist(cnf_locs):
            dip_mom_vec = cnf_fs[-1].file.dipole_moment.read(cnf_locs)
            dip_mom = automol.prop.total_dipole_moment(dip_mom_vec)
        else:
            dip_mom = None
        if cnf_fs[-1].file.polarizability.exist(cnf_locs):
            polar_tensor = cnf_fs[-1].file.polarizability.read(cnf_locs)
            polar = automol.prop.total_polarizability(polar_tensor)
        else:
            polar = None

        # Read the energy transfer filesystems
        if etrans_fs[-1].file.epsilon.exist(etrans_locs):
            eps = etrans_fs[-1].file.epsilon.read(etrans_locs)
        else:
            eps = None
        if etrans_fs[-1].file.sigma.exist(etrans_locs):
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
