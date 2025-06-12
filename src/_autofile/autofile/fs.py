""" This module generates managers for the various layers in the file system.

The overall file system layout is as follows::

    PREFIX/
        SPECIES/THEORY
            CONFORMER*/ [SINGLE POINT*/ or HIGH SPIN*/]
                SYMMETRY/
                ZMATRIX/
                    SCAN*/ [SINGLE POINT*/ or HIGH SPIN*/]
                    CSCAN*/ [SINGLE POINT*/ or HIGH SPIN*/]
        REACTION/THEORY/
            ZMATRIX/
                SCAN*/
            TRANSITION STATE*/
                CONFORMER*/ [SINGLE POINT*/ or HIGH SPIN*/]
                    SYMMETRY/
                    ZMATRIX/
                        SCAN*/ [SINGLE POINT*/ or HIGH SPIN*/]
                        CSCAN*/

    * Directory layer where electronic structure jobs are run.

Manager for the species-specific layer:
    - ``SPECIES``: :meth:`autofile.fs.species`

Managers for the reaction-specific layers:
    - ``REACTION``: :meth:`autofile.fs.reaction`
    - ``TRANSITION STATE``: :meth:`autofile.fs.transition_state`

Managers for layers used by both species and reaction file systems:
    - ``THEORY``: :meth:`autofile.fs.theory`
    - ``CONFORMER``: :meth:`autofile.fs.conformer`
    - ``SYMMETRY``: :meth:`autofile.fs.symmetry`
    - ``SINGLE POINT``: :meth:`autofile.fs.single_point`
    - ``HIGH SPIN``: :meth:`autofile.fs.high_spin`
    - ``ZMATRIX``: :meth:`autofile.fs.zmatrix`
    - ``SCAN``: :meth:`autofile.fs.scan`
    - ``CSCAN``: :meth:`autofile.fs.cscan`
    - ``VRCTST``: :meth:`autofile.fs.vrctst`

Each function in this module returns a tuple of autofile.model.DataSeries
objects for interacting with successive layers in a file system.
"""
import os
import pathlib
from autofile.schema import data_files
from autofile.schema import data_series
from autofile.schema import info_objects
from autofile.schema import json_objects


class _FilePrefix():
    """ file prefixes """
    RUN = 'run'
    BUILD = 'build'
    CONF = 'conf'
    RSAMP = 'rsamp'
    CSAMP = 'csamp'
    TAU = 'tau'
    SP = 'sp'
    HS = 'hs'
    SCAN = 'scan'
    GEOM = 'geom'
    GRAD = 'grad'
    HESS = 'hess'
    MIN = 'min'
    VPT2 = 'vpt2'
    FC = 'fc'
    LJ = 'lj'
    IRC = 'irc'
    ZMAT = 'zmat'
    INSTAB = 'instab'
    VRCTST = 'vrctst'


class FileAttributeName():
    """ DataFile attribute names """
    INFO = 'info'
    INFO2 = 'info2'
    INSTAB = 'instability'
    INPUT = 'input'
    OUTPUT = 'output'
    VMATRIX = 'vmatrix'
    TORS = 'torsions'
    RTORS = 'ring_torsions'
    GEOM_INFO = 'geometry_info'
    GRAD_INFO = 'gradient_info'
    HESS_INFO = 'hessian_info'
    VPT2_INFO = 'vpt2_info'
    IRC_INFO = 'irc_info'
    GEOM_INPUT = 'geometry_input'
    GRAD_INPUT = 'gradient_input'
    HESS_INPUT = 'hessian_input'
    VPT2_INPUT = 'vpt2_input'
    IRC_INPUT = 'irc_input'
    ENERGY = 'energy'
    INF_SEP_ENE = 'inf_sep_ene'
    GEOM = 'geometry'
    ZMAT = 'zmatrix'
    GRAD = 'gradient'
    HESS = 'hessian'
    HFREQ = 'harmonic_frequencies'
    TRAJ = 'trajectory'
    REAC = 'reaction'
    ANHFREQ = 'anharmonic_frequencies'
    ANHZPVE = 'anharmonic_zpve'
    CUBIC_FC = 'cubic_force_constants'
    QUARTIC_FC = 'quartic_force_constants'
    XMAT = 'anharmonicity_matrix'
    VIBROT_MAX = 'vibro_rot_alpha_matrix'
    CENTIF_DIST = 'quartic_centrifugal_dist_consts'
    LJ_EPS = 'lennard_jones_epsilon'
    LJ_SIG = 'lennard_jones_sigma'
    LJ_INPUT = 'lennard_jones_input'
    LJ_ELSTRUCT = 'lennard_jones_elstruct'
    LJ_INFO = 'lennard_jones_info'
    EXT_SYM = 'external_symmetry_number'
    INT_SYM = 'internal_symmetry_number'
    DIP_MOM = 'dipole_moment'
    POLAR = 'polarizability'
    VRC_TST = 'vrctst_tst'
    VRC_DIVSUR = 'vrctst_divsur'
    VRC_MOLP = 'vrctst_molpro'
    VRC_TML = 'vrctst_tml'
    VRC_STRUCT = 'vrctst_struct'
    VRC_POT = 'vrctst_pot'
    VRC_FLUX = 'vrctst_flux'


class _JSONAttributeName():
    """ JSON attribute names """
    INFO = 'info'
    INPUT = 'input'
    OUTPUT = 'output'
    VMATRIX = 'vmatrix'
    GEOM_INFO = 'geometry_info'
    GRAD_INFO = 'gradient_info'
    HESS_INFO = 'hessian_info'
    VPT2_INFO = 'vpt2_info'
    IRC_INFO = 'irc_info'
    GEOM_INPUT = 'geometry_input'
    GRAD_INPUT = 'gradient_input'
    HESS_INPUT = 'hessian_input'
    VPT2_INPUT = 'vpt2_input'
    IRC_INPUT = 'irc_input'
    ENERGY = 'energy'
    GEOM = 'geometry'
    ZMAT = 'zmatrix'
    GRAD = 'gradient'
    HESS = 'hessian'
    HFREQ = 'harmonic_frequencies'
    REAC = 'reaction'
    TRAJ = 'trajectory'
    ANHFREQ = 'anharmonic_frequencies'
    ANHZPVE = 'anharmonic_zpve'
    XMAT = 'anharmonicity_matrix'
    VIBROT_MAX = 'vibro_rot_alpha_matrix'
    CENTIF_DIST = 'quartic_centrifugal_dist_consts'
    LJ_EPS = 'lennard_jones_epsilon'
    LJ_SIG = 'lennard_jones_sigma'


class _LayerPrefix():
    """ file prefixes """
    RUN = 'RUN'
    SPC = 'SPC'
    CONF = 'CONF'
    RING = 'RING'
    TAU = 'TAU'
    SP = 'SP'


# Managers for species-specific layers
def species(prefix):
    """ Create a species file system manager for a given path.

    Layers:
        [0] The species "trunk" layer.
            Specifiers:
            []

            (no files)

            Generated by :meth:`autofile.schema.data_series.species_trunk`

        [1] The species "leaf" layer.
            Specifiers:
            [InChI String, Charge, Spin Multiplicity]

            (no files)

            Generated by :meth:`autofile.schema.data_series.species_leaf`


    :param prefix: Path to where the file system will be created
    :type prefix: str
    :returns: A tuple of DataSeries objects. Each item in this tuple manages a
        different layer in the file system, as described above.
    """
    trunk_ds = data_series.species_trunk(prefix)
    leaf_ds = data_series.species_leaf(prefix, root_ds=trunk_ds)
    return (trunk_ds, leaf_ds)


# Managers for reaction-specific layers
def reaction(prefix):
    """ construct the reaction filesystem (2 layers)

    locators:
        0 - []
                (no files)
        1 - [rxn_ichs, rxn_chgs, rxn_muls, ts_mul]
                (no files)

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = data_series.reaction_trunk(prefix)
    leaf_ds = data_series.reaction_leaf(prefix, root_ds=trunk_ds)

    return (trunk_ds, leaf_ds)


def transition_state(prefix):
    """ construct the ts filesystem (1 layer)

    locators:
        0 - []
                files:
                - energy
                - geometry
                - zmatrix
        1 - [TS Number]

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = data_series.transition_state_trunk(prefix)
    leaf_ds = data_series.transition_state_leaf(prefix, root_ds=trunk_ds)

    return (trunk_ds, leaf_ds)


# Managers for layers used by both species and reaction file systems
def theory(prefix):
    """ construct the theory filesystem (1 layer)

    locators:
        [0] "leaf" layer"
             Specifiers:
             [method, basis, orb_type]

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    leaf_ds = data_series.theory_leaf(prefix)
    return (leaf_ds,)


def conformer(prefix):
    """ Create a conformer file system manager for a given path.

    Layers:
        [0] The "trunk" layer.
            Specifiers:
            []

            Generated by :meth:`autofile.schema.data_series.conformer_trunk`

        [1] The "branch" layer.
            Specifiers:
            [Ring ID]

            Files:
                - info
                - energy
                - trajectory

            Generated by :meth:`autofile.schema.data_series.conformer_branch`

        [2] The "leaf" layer.
            Specifiers:
            [Conformer ID]

            Files:
                - geometry_info
                - gradient_info
                - hessian_info
                - geometry_input
                - gradient_input
                - hessian_input
                - geometry
                - gradient
                - hessian
                - harmonic_frequencies
                - vpt2_info
                - vpt2_input
                - anharmonic_frequencies
                - anharmonic_zpve
                - anharmonicity_matrix
                - vibro_rot_alpha_matrix
                - quartic_centrifugal_dist_consts
                - cubic_force_constants
                - quartic_force_constants
                - dipole_moment
                - polarizability

            Generated by :meth:`autofile.schema.data_series.conformer_leaf`


    :param prefix: Path to where the file system will be created
    :type prefix: str
    :returns: A tuple of DataSeries objects. Each item in this tuple manages a
        different layer in the file system, as described above.
    """
    trunk_ds = data_series.conformer_trunk(prefix)
    branch_ds = data_series.conformer_branch(prefix, root_ds=trunk_ds)
    leaf_ds = data_series.conformer_leaf(prefix, root_ds=branch_ds)

    # Trunk Files
    inf_dfile = data_files.information(
        _FilePrefix.RSAMP, function=info_objects.conformer_trunk)
    traj_dfile = data_files.trajectory(_FilePrefix.CONF)
    trunk_ds.add_data_files({
        FileAttributeName.INFO: inf_dfile,
        FileAttributeName.TRAJ: traj_dfile})

    # Branch Files
    inf2_dfile = data_files.information(
        _FilePrefix.CSAMP, function=info_objects.conformer_branch)
    traj2_dfile = data_files.trajectory(_FilePrefix.CONF)
    branch_ds.add_data_files({
        FileAttributeName.INFO: inf2_dfile,
        FileAttributeName.TRAJ: traj2_dfile})

    # Leaf files
    geom_inf_dfile = data_files.information(_FilePrefix.GEOM,
                                            function=info_objects.run)
    grad_inf_dfile = data_files.information(_FilePrefix.GRAD,
                                            function=info_objects.run)
    hess_inf_dfile = data_files.information(_FilePrefix.HESS,
                                            function=info_objects.run)
    vpt2_inf_dfile = data_files.information(_FilePrefix.VPT2,
                                            function=info_objects.vpt2)
    geom_inp_dfile = data_files.input_file(_FilePrefix.GEOM)
    grad_inp_dfile = data_files.input_file(_FilePrefix.GRAD)
    hess_inp_dfile = data_files.input_file(_FilePrefix.HESS)
    vpt2_inp_dfile = data_files.input_file(_FilePrefix.VPT2)
    geom_dfile = data_files.geometry(_FilePrefix.GEOM)
    grad_dfile = data_files.gradient(_FilePrefix.GRAD)
    hess_dfile = data_files.hessian(_FilePrefix.HESS)
    hfreq_dfile = data_files.harmonic_frequencies(_FilePrefix.HESS)
    anhfreq_dfile = data_files.anharmonic_frequencies(_FilePrefix.VPT2)
    anhzpve_dfile = data_files.anharmonic_zpve(_FilePrefix.VPT2)
    xmat_dfile = data_files.anharmonicity_matrix(_FilePrefix.VPT2)
    vibrot_mat_dfile = data_files.vibro_rot_alpha_matrix(_FilePrefix.VPT2)
    centrif_dist_dfile = data_files.quartic_centrifugal_dist_consts(
        _FilePrefix.VPT2)
    cubic_fc_dfile = data_files.cubic_force_constants(_FilePrefix.FC)
    quartic_fc_dfile = data_files.quartic_force_constants(_FilePrefix.FC)
    dip_mom_dfile = data_files.dipole_moment(_FilePrefix.GEOM)
    polar_dfile = data_files.polarizability(_FilePrefix.GEOM)

    leaf_ds.add_data_files({
        FileAttributeName.GEOM_INFO: geom_inf_dfile,
        FileAttributeName.GRAD_INFO: grad_inf_dfile,
        FileAttributeName.HESS_INFO: hess_inf_dfile,
        FileAttributeName.GEOM_INPUT: geom_inp_dfile,
        FileAttributeName.GRAD_INPUT: grad_inp_dfile,
        FileAttributeName.HESS_INPUT: hess_inp_dfile,
        FileAttributeName.GEOM: geom_dfile,
        FileAttributeName.GRAD: grad_dfile,
        FileAttributeName.HESS: hess_dfile,
        FileAttributeName.HFREQ: hfreq_dfile,
        FileAttributeName.VPT2_INFO: vpt2_inf_dfile,
        FileAttributeName.VPT2_INPUT: vpt2_inp_dfile,
        FileAttributeName.ANHFREQ: anhfreq_dfile,
        FileAttributeName.ANHZPVE: anhzpve_dfile,
        FileAttributeName.XMAT: xmat_dfile,
        FileAttributeName.VIBROT_MAX: vibrot_mat_dfile,
        FileAttributeName.CENTIF_DIST: centrif_dist_dfile,
        FileAttributeName.CUBIC_FC: cubic_fc_dfile,
        FileAttributeName.QUARTIC_FC: quartic_fc_dfile,
        FileAttributeName.DIP_MOM: dip_mom_dfile,
        FileAttributeName.POLAR: polar_dfile})

    return (trunk_ds, branch_ds, leaf_ds)


def single_point(prefix, json_layer=None):
    """ construct the single-point filesystem (2 layers)

    locators:
        0 - []
                (no files)
        1 - [method, basis, orb_type]
                files:
                - info
                - input
                - energy

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = data_series.single_point_trunk(prefix)
    leaf_ds = data_series.single_point_leaf(prefix, root_ds=trunk_ds)

    inp_dfile = data_files.input_file(_FilePrefix.SP)
    inf_dfile = data_files.information(_FilePrefix.SP,
                                       function=info_objects.run)
    ene_dfile = data_files.energy(_FilePrefix.SP)
    leaf_ds.add_data_files({
        FileAttributeName.INFO: inf_dfile,
        FileAttributeName.INPUT: inp_dfile,
        FileAttributeName.ENERGY: ene_dfile})
    inp_jobj = json_objects.input_file(
        _FilePrefix.SP, json_prefix=(json_layer, _LayerPrefix.SP))
    inf_jobj = json_objects.information(
        _FilePrefix.SP, json_prefix=(
            json_layer, _LayerPrefix.SP), function=info_objects.run)
    ene_jobj = json_objects.energy(
        _FilePrefix.SP, json_prefix=(json_layer, _LayerPrefix.SP))
    leaf_ds.add_json_entries({
        _JSONAttributeName.INFO: inf_jobj,
        _JSONAttributeName.INPUT: inp_jobj,
        _JSONAttributeName.ENERGY: ene_jobj})

    return (trunk_ds, leaf_ds)


def high_spin(prefix):
    """ construct the high-spin, single-point filesystem (2 layers)

    locators:
        0 - []
                (no files)
        1 - [method, basis, orb_type]
                files:
                - info
                - input
                - energy

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = data_series.high_spin_trunk(prefix)
    leaf_ds = data_series.high_spin_leaf(prefix, root_ds=trunk_ds)

    inp_dfile = data_files.input_file(_FilePrefix.HS)
    inf_dfile = data_files.information(_FilePrefix.HS,
                                       function=info_objects.run)
    ene_dfile = data_files.energy(_FilePrefix.HS)
    leaf_ds.add_data_files({
        FileAttributeName.INFO: inf_dfile,
        FileAttributeName.INPUT: inp_dfile,
        FileAttributeName.ENERGY: ene_dfile})

    return (trunk_ds, leaf_ds)


def symmetry(prefix):
    """ Create a symmetry-conformer file system manager for a given path.

    Layers:
        [0] The "trunk" layer.
            Specifiers:
            []

            Files:
                - trajectory

            Generated by :meth:`autofile.schema.data_series.symmetry_trunk`

        [1] The "leaf" layer.
            Specifiers:
            [Conformer ID]

            Files:
                - geometry

            Generated by :meth:`autofile.schema.data_series.conformer_leaf`


    :param prefix: Path to where the file system will be created
    :type prefix: str
    :returns: A tuple of DataSeries objects. Each item in this tuple manages a
        different layer in the file system, as described above.
    """
    trunk_ds = data_series.symmetry_trunk(prefix)
    leaf_ds = data_series.conformer_leaf(prefix, root_ds=trunk_ds)

    traj_dfile = data_files.trajectory(_FilePrefix.CONF)
    geom_dfile = data_files.geometry(_FilePrefix.GEOM)
    geom_inf_dfile = data_files.information(_FilePrefix.GEOM,
                                            function=info_objects.run)
    trunk_ds.add_data_files({
        FileAttributeName.TRAJ: traj_dfile})
    leaf_ds.add_data_files({
        FileAttributeName.GEOM: geom_dfile,
        FileAttributeName.GEOM_INFO: geom_inf_dfile})

    return (trunk_ds, leaf_ds)


def zmatrix(prefix):
    """ Create a z-matrix file system manager for a given path.

    Layers:
        [0] The zmatrix "trunk" layer.
            Specifiers:
            []

            (no files)

            Generated by :meth:`autofile.schema.data_series.zmatrix_trunk`

        [1] The zmatrix "leaf" layer.
            Specifiers:
            [Z-Matrix Number]

            Files:
                -

            Generated by :meth:`autofile.schema.data_series.zmatrix_leaf`


    :param prefix: Path to where the file system will be created
    :type prefix: str
    :returns: A tuple of DataSeries objects. Each item in this tuple manages a
        different layer in the file system, as described above.
    """
    trunk_ds = data_series.zmatrix_trunk(prefix)
    leaf_ds = data_series.zmatrix_leaf(prefix, root_ds=trunk_ds)

    instab_inf_dfile = data_files.instability(_FilePrefix.INSTAB)

    zmat_inf_dfile = data_files.information(
        _FilePrefix.ZMAT, function=info_objects.run)
    zmat_inp_dfile = data_files.input_file(_FilePrefix.ZMAT)

    tors_dfile = data_files.torsions(_FilePrefix.ZMAT)
    rtors_dfile = data_files.ring_torsions(_FilePrefix.ZMAT)

    zmat_dfile = data_files.zmatrix(_FilePrefix.ZMAT)
    reac_dfile = data_files.reaction(_FilePrefix.ZMAT)

    leaf_ds.add_data_files({
        FileAttributeName.INSTAB: instab_inf_dfile,
        FileAttributeName.GEOM_INFO:  zmat_inf_dfile,
        FileAttributeName.TORS: tors_dfile,
        FileAttributeName.RTORS: rtors_dfile,
        FileAttributeName.GEOM_INPUT: zmat_inp_dfile,
        FileAttributeName.ZMAT: zmat_dfile,
        FileAttributeName.REAC: reac_dfile})

    return (trunk_ds, leaf_ds)


def scan(prefix):
    """ construct the scan filesystem (3 layers)

    three layers with the following locators:
        0 - []
                files:
                - vmatrix
        1 - [coo_names]
                files:
                - info
                - trajectory
        2 - [coo_names, coo_vals]
                files:
                - geometry_info
                - gradient_info
                - hessian_info
                - irc_info
                - geometry_input
                - gradient_input
                - hessian_input
                - irc_input
                - energy
                - geometry
                - zmatrix
                - gradient
                - hessian
                - harmonic_frequencies

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = data_series.scan_trunk(prefix)
    branch_ds = data_series.scan_branch(prefix, root_ds=trunk_ds)
    leaf_ds = data_series.scan_leaf(prefix, root_ds=branch_ds)

    vma_dfile = data_files.vmatrix(_FilePrefix.SCAN)
    trunk_ds.add_data_files({
        FileAttributeName.VMATRIX: vma_dfile})

    inf_dfile = data_files.information(_FilePrefix.SCAN,
                                       function=info_objects.scan_branch)
    traj_dfile = data_files.trajectory(_FilePrefix.SCAN)
    branch_ds.add_data_files({
        FileAttributeName.INFO: inf_dfile,
        FileAttributeName.TRAJ: traj_dfile})

    # Need an irc file in the branch!
    # Need an run irc file in the forward and backward direction

    geom_inf_dfile = data_files.information(_FilePrefix.GEOM,
                                            function=info_objects.run)
    grad_inf_dfile = data_files.information(_FilePrefix.GRAD,
                                            function=info_objects.run)
    hess_inf_dfile = data_files.information(_FilePrefix.HESS,
                                            function=info_objects.run)
    irc_inf_dfile = data_files.information(_FilePrefix.IRC,
                                           function=info_objects.run)
    geom_inp_dfile = data_files.input_file(_FilePrefix.GEOM)
    grad_inp_dfile = data_files.input_file(_FilePrefix.GRAD)
    hess_inp_dfile = data_files.input_file(_FilePrefix.HESS)
    irc_inp_dfile = data_files.input_file(_FilePrefix.IRC)
    ene_dfile = data_files.energy(_FilePrefix.GEOM)
    geom_dfile = data_files.geometry(_FilePrefix.GEOM)
    zmat_dfile = data_files.zmatrix(_FilePrefix.GEOM)
    grad_dfile = data_files.gradient(_FilePrefix.GRAD)
    hess_dfile = data_files.hessian(_FilePrefix.HESS)
    hfreq_dfile = data_files.harmonic_frequencies(_FilePrefix.HESS)
    leaf_ds.add_data_files({
        FileAttributeName.GEOM_INFO: geom_inf_dfile,
        FileAttributeName.GRAD_INFO: grad_inf_dfile,
        FileAttributeName.HESS_INFO: hess_inf_dfile,
        FileAttributeName.IRC_INFO: irc_inf_dfile,
        FileAttributeName.GEOM_INPUT: geom_inp_dfile,
        FileAttributeName.GRAD_INPUT: grad_inp_dfile,
        FileAttributeName.HESS_INPUT: hess_inp_dfile,
        FileAttributeName.IRC_INPUT: irc_inp_dfile,
        FileAttributeName.ENERGY: ene_dfile,
        FileAttributeName.GEOM: geom_dfile,
        FileAttributeName.ZMAT: zmat_dfile,
        FileAttributeName.GRAD: grad_dfile,
        FileAttributeName.HESS: hess_dfile,
        FileAttributeName.HFREQ: hfreq_dfile})

    return (trunk_ds, branch_ds, leaf_ds)


def cscan(prefix):
    """ Create a constrained scan (c-scan) file system manager for a given path.

    Layers:
        [0] The c-scan "trunk" layer.
            Specifiers:
            []

            (no files)

            Generated by :meth:`autofile.schema.data_series.cscan_trunk`

        [1] The first c-scan "branch" layer.
            Specifiers:
            [Constraints]

            Files:
                - info
                - trajectory

            Generated by :meth:`autofile.schema.data_series.cscan_branch1`

        [2] The second c-scan "branch" layer.
            Specifiers:
            [Constraints, Scan Coordinate Names]

            (no files)

            Generated by :meth:`autofile.schema.data_series.cscan_branch2`

        [3] The c-scan "leaf" layer.
            Specifiers:
            [Constraints, Scan Coordinate Names, Scan Coordinate Values]

            Files:
                - geometry_info
                - gradient_info
                - hessian_info
                - geometry_input
                - gradient_input
                - hessian_input
                - energy
                - geometry
                - zmatrix
                - gradient
                - hessian
                - harmonic_frequencies

            Generated by :meth:`autofile.schema.data_series.cscan_leaf`


    :param prefix: Path to where the file system will be created
    :type prefix: str
    :returns: A tuple of DataSeries objects. Each item in this tuple manages a
        different layer in the file system, as described above.
    """
    trunk_ds = data_series.cscan_trunk(prefix)
    branch1_ds = data_series.cscan_branch1(prefix, root_ds=trunk_ds)
    branch2_ds = data_series.cscan_branch2(prefix, root_ds=branch1_ds)
    leaf_ds = data_series.cscan_leaf(prefix, root_ds=branch2_ds)

    inf_dfile = data_files.information(_FilePrefix.SCAN,
                                       function=info_objects.scan_branch)
    traj_dfile = data_files.trajectory(_FilePrefix.SCAN)

    branch1_ds.add_data_files({
        FileAttributeName.INFO: inf_dfile,
        FileAttributeName.TRAJ: traj_dfile})

    inf_sep_ene_dfile = data_files.energy(_FilePrefix.VRCTST)

    branch2_ds.add_data_files({
        FileAttributeName.INF_SEP_ENE: inf_sep_ene_dfile})

    geom_inf_dfile = data_files.information(_FilePrefix.GEOM,
                                            function=info_objects.run)
    grad_inf_dfile = data_files.information(_FilePrefix.GRAD,
                                            function=info_objects.run)
    hess_inf_dfile = data_files.information(_FilePrefix.HESS,
                                            function=info_objects.run)
    geom_inp_dfile = data_files.input_file(_FilePrefix.GEOM)
    grad_inp_dfile = data_files.input_file(_FilePrefix.GRAD)
    hess_inp_dfile = data_files.input_file(_FilePrefix.HESS)
    ene_dfile = data_files.energy(_FilePrefix.GEOM)
    geom_dfile = data_files.geometry(_FilePrefix.GEOM)
    zmat_dfile = data_files.zmatrix(_FilePrefix.GEOM)
    grad_dfile = data_files.gradient(_FilePrefix.GRAD)
    hess_dfile = data_files.hessian(_FilePrefix.HESS)
    hfreq_dfile = data_files.harmonic_frequencies(_FilePrefix.HESS)
    leaf_ds.add_data_files({
        FileAttributeName.GEOM_INFO: geom_inf_dfile,
        FileAttributeName.GRAD_INFO: grad_inf_dfile,
        FileAttributeName.HESS_INFO: hess_inf_dfile,
        FileAttributeName.GEOM_INPUT: geom_inp_dfile,
        FileAttributeName.GRAD_INPUT: grad_inp_dfile,
        FileAttributeName.HESS_INPUT: hess_inp_dfile,
        FileAttributeName.ENERGY: ene_dfile,
        FileAttributeName.GEOM: geom_dfile,
        FileAttributeName.ZMAT: zmat_dfile,
        FileAttributeName.GRAD: grad_dfile,
        FileAttributeName.HESS: hess_dfile,
        FileAttributeName.HFREQ: hfreq_dfile})

    return (trunk_ds, branch1_ds, branch2_ds, leaf_ds)


def tau(prefix):
    """ construct the tau filesystem (2 layers)

    locators:
        0 - []
                files:
                - vmatrix
                - info
                - trajectory
        0 - [tid]
                files:
                - geometry_info
                - gradient_info
                - hessian_info
                - geometry_input
                - gradient_input
                - hessian_input
                - energy
                - geometry
                - zmatrix
                - gradient
                - hessian
                - harmonic_frequencies

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = data_series.tau_trunk(prefix)
    leaf_ds = data_series.tau_leaf(prefix, root_ds=trunk_ds)

    vma_dfile = data_files.vmatrix(_FilePrefix.TAU)
    inf_dfile = data_files.information(_FilePrefix.TAU,
                                       function=info_objects.tau_trunk)
    traj_dfile = data_files.trajectory(_FilePrefix.TAU)
    trunk_ds.add_data_files({
        FileAttributeName.VMATRIX: vma_dfile,
        FileAttributeName.INFO: inf_dfile,
        FileAttributeName.TRAJ: traj_dfile})

    geom_inf_dfile = data_files.information(_FilePrefix.GEOM,
                                            function=info_objects.run)
    grad_inf_dfile = data_files.information(_FilePrefix.GRAD,
                                            function=info_objects.run)
    hess_inf_dfile = data_files.information(_FilePrefix.HESS,
                                            function=info_objects.run)
    geom_inp_dfile = data_files.input_file(_FilePrefix.GEOM)
    grad_inp_dfile = data_files.input_file(_FilePrefix.GRAD)
    hess_inp_dfile = data_files.input_file(_FilePrefix.HESS)
    ene_dfile = data_files.energy(_FilePrefix.GEOM)
    geom_dfile = data_files.geometry(_FilePrefix.GEOM)
    zmat_dfile = data_files.zmatrix(_FilePrefix.GEOM)
    grad_dfile = data_files.gradient(_FilePrefix.GRAD)
    hess_dfile = data_files.hessian(_FilePrefix.HESS)
    hfreq_dfile = data_files.harmonic_frequencies(_FilePrefix.HESS)
    leaf_ds.add_data_files({
        FileAttributeName.GEOM_INFO: geom_inf_dfile,
        FileAttributeName.GRAD_INFO: grad_inf_dfile,
        FileAttributeName.HESS_INFO: hess_inf_dfile,
        FileAttributeName.GEOM_INPUT: geom_inp_dfile,
        FileAttributeName.GRAD_INPUT: grad_inp_dfile,
        FileAttributeName.HESS_INPUT: hess_inp_dfile,
        FileAttributeName.ENERGY: ene_dfile,
        FileAttributeName.GEOM: geom_dfile,
        FileAttributeName.ZMAT: zmat_dfile,
        FileAttributeName.GRAD: grad_dfile,
        FileAttributeName.HESS: hess_dfile,
        FileAttributeName.HFREQ: hfreq_dfile})

    geom_inf_jobj = json_objects.information(
        _FilePrefix.GEOM, function=info_objects.run)
    grad_inf_jobj = json_objects.information(
        _FilePrefix.GRAD, function=info_objects.run)
    hess_inf_jobj = json_objects.information(
        _FilePrefix.HESS, function=info_objects.run)
    geom_inp_jobj = json_objects.input_file(_FilePrefix.GEOM)
    grad_inp_jobj = json_objects.input_file(_FilePrefix.GRAD)
    hess_inp_jobj = json_objects.input_file(_FilePrefix.HESS)
    ene_jobj = json_objects.energy(_FilePrefix.GEOM)
    geom_jobj = json_objects.geometry(_FilePrefix.GEOM)
    zmat_jobj = json_objects.zmatrix(_FilePrefix.GEOM)
    grad_jobj = json_objects.gradient(_FilePrefix.GRAD)
    hess_jobj = json_objects.hessian(_FilePrefix.HESS)
    hfreq_jobj = json_objects.harmonic_frequencies(_FilePrefix.HESS)

    leaf_ds.add_json_entries({
        _JSONAttributeName.GEOM_INFO: geom_inf_jobj,
        _JSONAttributeName.GRAD_INFO: grad_inf_jobj,
        _JSONAttributeName.HESS_INFO: hess_inf_jobj,
        _JSONAttributeName.GEOM_INPUT: geom_inp_jobj,
        _JSONAttributeName.GRAD_INPUT: grad_inp_jobj,
        _JSONAttributeName.HESS_INPUT: hess_inp_jobj,
        _JSONAttributeName.ENERGY: ene_jobj,
        _JSONAttributeName.GEOM: geom_jobj,
        _JSONAttributeName.ZMAT: zmat_jobj,
        _JSONAttributeName.GRAD: grad_jobj,
        _JSONAttributeName.HESS: hess_jobj,
        _JSONAttributeName.HFREQ: hfreq_jobj})

    return (trunk_ds, leaf_ds)


def energy_transfer(prefix):
    """ construct the energy transfer filesystem (3 layers)

    locators:
        0 - []
                files:
                - info
        1 - [ich, chg, mul]
                (no files)
        2 - [ich, chg, mul, method, basis, orb_type]
                files:
                - lennard_jones_epsilon
                - lennard_jones_sigma
                - trajectory

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = data_series.energy_transfer_trunk(prefix)
    branch_ds = data_series.energy_transfer_branch(prefix, root_ds=trunk_ds)
    leaf_ds = data_series.energy_transfer_leaf(prefix, root_ds=branch_ds)

    inf_dfile = data_files.information(
        _FilePrefix.LJ, function=info_objects.lennard_jones)
    lj_inp_dfile = data_files.lennard_jones_input(_FilePrefix.LJ)
    temp_inp_dfile = data_files.lennard_jones_elstruct(_FilePrefix.LJ)

    eps_dfile = data_files.lennard_jones_epsilon(_FilePrefix.LJ)
    sig_dfile = data_files.lennard_jones_sigma(_FilePrefix.LJ)
    traj_dfile = data_files.trajectory(_FilePrefix.LJ)

    leaf_ds.add_data_files({
        FileAttributeName.LJ_INFO: inf_dfile,
        FileAttributeName.LJ_INPUT: lj_inp_dfile,
        FileAttributeName.LJ_ELSTRUCT: temp_inp_dfile,
        FileAttributeName.LJ_EPS: eps_dfile,
        FileAttributeName.LJ_SIG: sig_dfile,
        FileAttributeName.TRAJ: traj_dfile})

    return (trunk_ds, branch_ds, leaf_ds)


def vrctst(prefix):
    """ Create a VRC-TST file system manager for a given path.

    Layers:
        [0] The VRC-TST "trunk" layer.
            Specifiers:
            []

            (no files)

            Generated by :meth:`autofile.schema.data_series.vrctst_trunk`

        [1] The VRC-TST "leaf" layer.
            Specifiers:
            [VRC-TST Number]

            Files:
                - tst.inp (VaReCoF Input)
                - divsur.inp (VaReCoF Input)
                - els.tml (VaReCoF Input)
                - structure.inp (VaReCoF Input)
                - corr_pot.f (VaReCoF Input)
                - flux file

            Generated by :meth:`autofile.schema.data_series.vrctst_leaf`


    :param prefix: Path to where the file system will be created
    :type prefix: str
    :returns: A tuple of DataSeries objects. Each item in this tuple manages a
        different layer in the file system, as described above.
    """

    trunk_ds = data_series.vrctst_trunk(prefix)
    leaf_ds = data_series.vrctst_leaf(prefix, root_ds=trunk_ds)

    tst_inp_dfile = data_files.vrctst_tst(_FilePrefix.VRCTST)
    divsur_inp_dfile = data_files.vrctst_divsur(_FilePrefix.VRCTST)
    molp_inp_dfile = data_files.vrctst_molpro(_FilePrefix.VRCTST)
    tml_inp_dfile = data_files.vrctst_tml(_FilePrefix.VRCTST)
    struct_inp_dfile = data_files.vrctst_struct(_FilePrefix.VRCTST)
    pot_inp_dfile = data_files.vrctst_pot(_FilePrefix.VRCTST)
    flux_dfile = data_files.vrctst_flux(_FilePrefix.VRCTST)

    leaf_ds.add_data_files({
        FileAttributeName.VRC_TST: tst_inp_dfile,
        FileAttributeName.VRC_DIVSUR: divsur_inp_dfile,
        FileAttributeName.VRC_MOLP: molp_inp_dfile,
        FileAttributeName.VRC_TML: tml_inp_dfile,
        FileAttributeName.VRC_STRUCT: struct_inp_dfile,
        FileAttributeName.VRC_POT: pot_inp_dfile,
        FileAttributeName.VRC_FLUX: flux_dfile
    })

    return (trunk_ds, leaf_ds)


# Managers specific to the run file system
def run(prefix):
    """ construct the run filesystem (2 layers)

    locators:
        0 - []
                files:
                - info
        1 - [job]
                files:
                - info
                - input
                - output
                - geometry
                - zmatrix

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = data_series.run_trunk(prefix)
    leaf_ds = data_series.run_leaf(prefix, root_ds=trunk_ds)

    inf_dfile = data_files.information(_FilePrefix.RUN,
                                       function=info_objects.run)
    inp_dfile = data_files.input_file(_FilePrefix.RUN)
    out_dfile = data_files.output_file(_FilePrefix.RUN)
    geom_dfile = data_files.geometry(_FilePrefix.GEOM)
    zmat_dfile = data_files.zmatrix(_FilePrefix.ZMAT)

    trunk_ds.add_data_files({
        FileAttributeName.INFO: inf_dfile})
    leaf_ds.add_data_files({
        FileAttributeName.INFO: inf_dfile,
        FileAttributeName.INPUT: inp_dfile,
        FileAttributeName.OUTPUT: out_dfile,
        FileAttributeName.GEOM: geom_dfile,
        FileAttributeName.ZMAT: zmat_dfile})

    return (trunk_ds, leaf_ds)


def subrun(prefix):
    """ construct the subrun filesystem (1 layer)

    locators:
        0 - [macro_idx, micro_idx]
                files:
                - info
                - input
                - output

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    leaf_ds = data_series.subrun_leaf(prefix)

    inf_dfile = data_files.information(_FilePrefix.RUN,
                                       function=info_objects.run)
    inp_dfile = data_files.input_file(_FilePrefix.RUN)
    out_dfile = data_files.output_file(_FilePrefix.RUN)
    leaf_ds.add_data_files({
        FileAttributeName.INFO: inf_dfile,
        FileAttributeName.INPUT: inp_dfile,
        FileAttributeName.OUTPUT: out_dfile})

    return (leaf_ds,)


def build(prefix):
    """ construct the build filesystem (2 layers)

    locators:
        0 - [job]
                (no files)
        2 - [job, formula]
                (no files)
        2 - [job, formula, num]
                files:
                - input
                - output

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = data_series.build_trunk(prefix)
    branch_ds = data_series.build_branch(prefix, root_ds=trunk_ds)
    leaf_ds = data_series.build_leaf(prefix, root_ds=branch_ds)

    inp_dfile = data_files.input_file(_FilePrefix.BUILD)
    out_dfile = data_files.output_file(_FilePrefix.BUILD)
    leaf_ds.add_data_files({
        FileAttributeName.INPUT: inp_dfile,
        FileAttributeName.OUTPUT: out_dfile})

    return (trunk_ds, branch_ds, leaf_ds)


def _process_root_args(root_fs=None, top_ds_name=None):
    if root_fs is not None:
        root_fs = dict(root_fs)
        assert top_ds_name in root_fs
        top_dsdir = root_fs[top_ds_name].dir
    else:
        root_fs = {}
        top_dsdir = None
    return root_fs, top_dsdir


FILE_SYSTEM_MANAGER_DCT = {
    'SPECIES': species,
    'REACTION': reaction,
    'TRANSITION STATE': transition_state,
    'THEORY': theory,
    'CONFORMER': conformer,
    'SYMMETRY': symmetry,
    'SINGLE POINT': single_point,
    'HIGH SPIN': high_spin,
    'ZMATRIX': zmatrix,
    'SCAN': scan,
    'CSCAN': cscan,
    'TAU': tau,
    'ENERGY TRANSFER': energy_transfer,
    'VRCTST': vrctst,
    'RUN': run,
    'SRUN': subrun,
    'BUILD': build
}


def path(pfx, key_locs_lst):
    """ Get the path through a file system hierarchy
    """
    pth = pfx
    for key_locs in key_locs_lst:
        assert len(key_locs) == 2
        key, locs = key_locs

        assert key in FILE_SYSTEM_MANAGER_DCT

        fs_ = FILE_SYSTEM_MANAGER_DCT[key](pth)
        fs_[-1].create(locs)  # run create command to filesys fix
        pth = os.path.join(pth, fs_[-1].path(locs))

    return pth


def manager(pfx, key_locs_lst, key):
    """ Get the manager for a specific part of the file system
    """
    pth = path(pfx, key_locs_lst)
    fs_ = _manager(pth, key)
    return fs_


def _manager(pfx, key):
    """ Get the manager for a specific part of the file system
    """
    fs_ = FILE_SYSTEM_MANAGER_DCT[key](pfx)
    return fs_


def iterate_locators(pfx, keys):
    """ Iterate over locators for all existing paths
    """
    depth = len(keys)
    locs_lst = [None] * depth

    def _iterate_locators(pfx_, keys_):
        if len(keys_) == 1:
            key_, = keys_

            fs_ = _manager(pfx_, key_)
            for locs in fs_[-1].existing():
                locs_lst[-1] = locs
                yield tuple(locs_lst)
        else:
            idx = depth - len(keys_)
            key_, keys_ = keys_[0], keys_[1:]

            fs_ = _manager(pfx_, key_)
            for locs in fs_[-1].existing():
                pfx_ = fs_[-1].path(locs)
                locs_lst[idx] = locs
                yield from _iterate_locators(pfx_, keys_)

    yield from _iterate_locators(pfx, keys)


def iterate_paths(pfx, keys):
    """ Iterate over all existing paths
    """
    if len(keys) == 1:
        key, = keys

        fs_ = _manager(pfx, key)
        for locs in fs_[-1].existing():
            yield fs_[-1].path(locs)
    else:
        key, keys = keys[0], keys[1:]
        fs_ = _manager(pfx, key)
        for locs in fs_[-1].existing():
            pfx_ = fs_[-1].path(locs)
            yield from iterate_paths(pfx_, keys)


def iterate_managers(pfx, keys, key):
    """ Iterate over managers at a specific level in the file system hierarchy
    """
    for pth in iterate_paths(pfx, keys):
        yield _manager(pth, key)


# Extra path manipulations
def path_prefix(pth, keys):
    """ Given a path and some layer keys, find the prefix

        :param pth: The path to some point in the file system
        :param keys: Keys to layers on top of the desired prefix
        :returns: The prefix path (the original path, without the layers
            specified by keys)
    """
    mgrs = [_manager('', key) for key in keys]
    depth = sum(ds.depth for m in mgrs for ds in m)
    pfx = os.path.join(*pathlib.PurePath(pth).parts[:-depth])
    return pfx
