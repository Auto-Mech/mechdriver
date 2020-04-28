"""
tau stuff
"""

import mess_io
import automol
import autofile
from routines.pf.messf import _tors as tors
from lib import filesys
from lib.phydat import phycon


def write_monte_carlo_mess_strings(tors_min_cnf_locs, tors_cnf_save_fs,
                                   spc_dct_i,
                                   frm_bnd_key, brk_bnd_key,
                                   elec_levels,
                                   tau_dat_file_name,
                                   freqs=(),
                                   saddle=False):
    """ Write out the input string for tau samling
    """
    # Loop over the torsions to get the flux strings
    flux_mode_str = ''
    if tors_min_cnf_locs is not None:

        # Get geometry for the torsional minimum
        zma = tors_cnf_save_fs[-1].file.zmatrix.read(
            tors_min_cnf_locs)
        name_matrix = automol.zmatrix.name_matrix(zma)
        key_matrix = automol.zmatrix.key_matrix(zma)
        tors_geo = tors_cnf_save_fs[-1].file.geometry.read(
            tors_min_cnf_locs)

        # Set torsional stuff
        tors_names = tors.get_tors_names(
            spc_dct_i, tors_cnf_save_fs, saddle=saddle)
        tors_sym_nums = tors.get_tors_sym_nums(
            spc_dct_i, tors_min_cnf_locs, tors_cnf_save_fs,
            frm_bnd_key, brk_bnd_key, saddle=False)

        # Write the MESS flux strings for each of the modes
        for tors_name, tors_sym in zip(tors_names, tors_sym_nums):

            # Get the idxs for the atoms used to define the dihedrals
            # Move code at some point to automol
            mode_idxs = [0, 0, 0, 0]
            for i, name_row in enumerate(name_matrix):
                if tors_name in name_row:
                    mode_idxs[0] = i
                    mode_idxs[1], mode_idxs[2], mode_idxs[3] = key_matrix[i]
                    break
            mode_idxs = [idx+1 for idx in mode_idxs]

            # Determine the flux span using the symmetry number
            mode_span = 360.0 / tors_sym

            # Write flux string for mode_ to overall flux string
            flux_mode_str += mess_io.writer.fluxional_mode(
                mode_idxs, span=mode_span)

    # Write the core string (seperate energies?)
    ground_ene = -0.02
    reference_ene = 0.00
    monte_carlo_str = mess_io.writer.monte_carlo(
        tors_geo, elec_levels,
        flux_mode_str, tau_dat_file_name,
        ground_ene, reference_ene,
        freqs=freqs, no_qc_corr=True, use_cm_shift=True)

    return monte_carlo_str


def write_tau_data_str(
        name, save_prefix,
        gradient=False, hessian=False):
    """ Write out data fle for partition function evaluation
    """
    cnf_save_fs = autofile.fs.conformer(save_prefix)
    min_cnf_locs = filesys.mincnfnf.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        ene_ref = cnf_save_fs[-1].file.energy.read(min_cnf_locs)

    tau_save_fs = autofile.fs.tau(save_prefix)
    evr = name+'\n'
    # cycle through saved tau geometries
    idx = 0
    for locs in tau_save_fs[-1].existing():
        geo = tau_save_fs[-1].file.geometry.read(locs)
        ene = tau_save_fs[-1].file.energy.read(locs)
        ene = (ene - ene_ref) * phycon.EH2KCAL
        ene_str = autofile.file.write.energy(ene)
        geo_str = autofile.file.write.geometry(geo)

        idx += 1
        idx_str = str(idx)

        evr += 'Sampling point'+idx_str+'\n'
        evr += 'Energy'+'\n'
        evr += ene_str+'\n'
        evr += 'Geometry'+'\n'
        evr += geo_str+'\n'
        if gradient:
            grad = tau_save_fs[-1].file.gradient.read(locs)
            grad_str = autofile.file.write.gradient(grad)
            evr += 'Gradient'+'\n'
            evr += grad_str
        if hessian:
            hess = tau_save_fs[-1].file.hessian.read(locs)
            hess_str = autofile.file.write.hessian(hess)
            evr += 'Hessian'+'\n'
            evr += hess_str+'\n'

    return evr
