"""
  Tests the varecof_io.writer functions
"""

import varecof_io.writer

# required arguments
mep_distances = [1.5958, 1.6958, 1.7958, 1.8958, 1.9958, 2.0958, 2.1958, 2.2958, 2.3958, 2.4958]
potentials = [
    [0.052, 0.175, 0.430, 0.724, 0.996, 1.199, 1.308, 1.317, 1.243, 1.113],
    [-0.722, -0.517, -0.372, -0.277, -0.218, -0.181, -0.153, -0.126, -0.096, -0.064],
    [0.224, 0.329, 0.556, 0.823, 1.071, 1.255, 1.348, 1.346, 1.263, 1.127]
]
bond_frm_idxs = [1, 3]
fortran_compiler = 'gfortran'
fortran_mak_path = '.'

# optional arguments
bnd_frm_syms = ['C', 'O']
species_name = 'mol'
pot_labels = ['basis+relaxed', 'basis', 'relaxed']
pot_file_names = ['mol']


def build_correction_potential(
    mep_vals,
    potentials,
    bond_form_idxs,
    fortran_compiler,
    fortran_make_path,
    bond_form_syms=(),
    species_name='',
    pot_labels=(),
    pot_file_names=()):
    """  use the MEP potentials to compile the correction potential .so file
    """

    # Write strings corresponding to each of the correction potential files
    species_corr_str = varecof_io.writer.corr_potentials.species(
        mep_distances,
        potentials,
        bnd_frm_idxs,
        bnd_frm_syms=bnd_frm_syms,
        species_name=species_name, 
        pot_labels=pot_labels)
    dummy_corr_str = varecof_io.writer.corr_potentials.dummy()
    pot_aux_str = varecof_io.writer.corr_potentials.auxiliary()
    makefile_str = varecof_io.writer.corr_potentials.makefile(
        fortran_compiler,
        pot_file_names=pot_file_names)

    # Write all of the files needed to build the correction potential
    with open('mol_corr.f', 'w') as mol_corr_file:
        mol_corr_file.write(species_corr_str)
    with open('dummy_corr.f', 'w') as dummy_corr_file:
        dummy_corr_file.write(dummy_corr_str)
    with open('pot_aux.f', 'w') as pot_aux_file:
        pot_aux_file.write(pot_aux_str)
    with open('makefile', 'w') as makefile_file:
        makefile_file.write(makefile_str)

    # Compile the correction potential
    varecof_io.writer.corr_potentials.compile_corr_pot(
        fortran_make_path)
