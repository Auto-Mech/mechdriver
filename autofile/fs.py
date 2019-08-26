""" filesystem library
"""
from autofile.system import file_
from autofile.system import info
from autofile.system import dir_
from autofile.system import model


class FilePrefix():
    """ file prefixes """
    RUN = 'run'
    BUILD = 'build'
    CONF = 'conf'
    TAU = 'tau'
    SP = 'sp'
    SCAN = 'scan'
    GEOM = 'geom'
    GRAD = 'grad'
    HESS = 'hess'
    MIN = 'min'
    VPT2 = 'vpt2'


class FileAttributeName():
    """ DataFile attribute names """
    INFO = 'info'
    INPUT = 'input'
    OUTPUT = 'output'
    VMATRIX = 'vmatrix'
    GEOM_INFO = 'geometry_info'
    GRAD_INFO = 'gradient_info'
    HESS_INFO = 'hessian_info'
    VPT2_INFO = 'vpt2_info'
    GEOM_INPUT = 'geometry_input'
    GRAD_INPUT = 'gradient_input'
    HESS_INPUT = 'hessian_input'
    VPT2_INPUT = 'vpt2_input'
    ENERGY = 'energy'
    GEOM = 'geometry'
    ZMAT = 'zmatrix'
    GRAD = 'gradient'
    HESS = 'hessian'
    HFREQ = 'harmonic_frequencies'
    TRAJ = 'trajectory'
    XMAT = 'anharmonicity_matrix'


class SeriesAttributeName():
    """ DataSeries attribute names
    """
    TRUNK = 'trunk'
    BRANCH = 'branch'
    LEAF = 'leaf'


def species(prefix):
    """ construct the species filesystem [trunk/leaf]

    layers:
     - trunk (specifiers: [])
     - leaf (specifiers: [ich, chg, mul])

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = dir_.species_trunk(prefix)
    leaf_ds = dir_.species_leaf(prefix, root_ds=trunk_ds)

    dir_fs = model.FileSystem({SeriesAttributeName.TRUNK: trunk_ds,
                               SeriesAttributeName.LEAF: leaf_ds})
    return dir_fs


def theory(prefix):
    """ construct the theory filesystem [leaf]

    layers:
     - leaf (specifiers: [method, basis, orb_restricted])

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    leaf_ds = dir_.theory_leaf(prefix)

    geom_dfile = file_.geometry(FilePrefix.GEOM)
    ene_dfile = file_.energy(FilePrefix.GEOM)
    zmat_dfile = file_.zmatrix(FilePrefix.GEOM)
    hess_dfile = file_.hessian(FilePrefix.HESS)
    leaf_ds.add_data_files({
        FileAttributeName.ENERGY: ene_dfile,
        FileAttributeName.GEOM: geom_dfile,
        FileAttributeName.HESS: hess_dfile,
        FileAttributeName.ZMAT: zmat_dfile})

    dir_fs = model.FileSystem({SeriesAttributeName.LEAF: leaf_ds})
    return dir_fs


def conformer(prefix):
    """ construct the conformer filesystem [trunk/leaf]

    layers:
     - trunk (specifiers: [])
     - leaf (specifiers: [cid])

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = dir_.conformer_trunk(prefix)
    leaf_ds = dir_.conformer_leaf(prefix, root_ds=trunk_ds)

    min_ene_dfile = file_.energy(FilePrefix.MIN)
    vma_dfile = file_.vmatrix(FilePrefix.CONF)
    inf_dfile = file_.information(FilePrefix.CONF,
                                  function=info.conformer_trunk)
    traj_dfile = file_.trajectory(FilePrefix.CONF)
    trunk_ds.add_data_files({
        FileAttributeName.VMATRIX: vma_dfile,
        FileAttributeName.INFO: inf_dfile,
        FileAttributeName.ENERGY: min_ene_dfile,
        FileAttributeName.TRAJ: traj_dfile})

    geom_inf_dfile = file_.information(FilePrefix.GEOM, function=info.run)
    grad_inf_dfile = file_.information(FilePrefix.GRAD, function=info.run)
    hess_inf_dfile = file_.information(FilePrefix.HESS, function=info.run)
    geom_inp_dfile = file_.input_file(FilePrefix.GEOM)
    grad_inp_dfile = file_.input_file(FilePrefix.GRAD)
    hess_inp_dfile = file_.input_file(FilePrefix.HESS)
    ene_dfile = file_.energy(FilePrefix.GEOM)
    geom_dfile = file_.geometry(FilePrefix.GEOM)
    zmat_dfile = file_.zmatrix(FilePrefix.GEOM)
    grad_dfile = file_.gradient(FilePrefix.GRAD)
    hess_dfile = file_.hessian(FilePrefix.HESS)
    hfreq_dfile = file_.harmonic_frequencies(FilePrefix.HESS)
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

    dir_fs = model.FileSystem({SeriesAttributeName.TRUNK: trunk_ds,
                               SeriesAttributeName.LEAF: leaf_ds})
    return dir_fs


def single_point(prefix):
    """ construct the single-point filesystem [trunk/leaf]

    layers:
     - trunk (specifiers: [])
     - leaf (specifiers: [method, basis, orb_restricted])

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = dir_.single_point_trunk(prefix)
    leaf_ds = dir_.single_point_leaf(prefix, root_ds=trunk_ds)

    inp_dfile = file_.input_file(FilePrefix.SP)
    inf_dfile = file_.information(FilePrefix.SP, function=info.run)
    ene_dfile = file_.energy(FilePrefix.SP)
    leaf_ds.add_data_files({
        FileAttributeName.INFO: inf_dfile,
        FileAttributeName.INPUT: inp_dfile,
        FileAttributeName.ENERGY: ene_dfile})

    dir_fs = model.FileSystem({
        SeriesAttributeName.TRUNK: trunk_ds,
        SeriesAttributeName.LEAF: leaf_ds})
    return dir_fs


def scan(prefix):
    """ construct the scan filesystem

    layers:
     - trunk (specifiers: [])
     - branch (specifiers: [coo_names])
     - leaf (specifiers: [coo_names, grid_idxs])

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = dir_.scan_trunk(prefix)
    branch_ds = dir_.scan_branch(prefix, root_ds=trunk_ds)
    leaf_ds = dir_.scan_leaf(prefix, root_ds=branch_ds)

    vma_dfile = file_.vmatrix(FilePrefix.SCAN)
    trunk_ds.add_data_files({
        FileAttributeName.VMATRIX: vma_dfile})

    inf_dfile = file_.information(FilePrefix.SCAN, function=info.scan_branch)
    traj_dfile = file_.trajectory(FilePrefix.SCAN)
    branch_ds.add_data_files({
        FileAttributeName.INFO: inf_dfile,
        FileAttributeName.TRAJ: traj_dfile})

    geom_inf_dfile = file_.information(FilePrefix.GEOM, function=info.run)
    grad_inf_dfile = file_.information(FilePrefix.GRAD, function=info.run)
    hess_inf_dfile = file_.information(FilePrefix.HESS, function=info.run)
    geom_inp_dfile = file_.input_file(FilePrefix.GEOM)
    grad_inp_dfile = file_.input_file(FilePrefix.GRAD)
    hess_inp_dfile = file_.input_file(FilePrefix.HESS)
    ene_dfile = file_.energy(FilePrefix.GEOM)
    geom_dfile = file_.geometry(FilePrefix.GEOM)
    zmat_dfile = file_.zmatrix(FilePrefix.GEOM)
    grad_dfile = file_.gradient(FilePrefix.GRAD)
    hess_dfile = file_.hessian(FilePrefix.HESS)
    hfreq_dfile = file_.harmonic_frequencies(FilePrefix.HESS)
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

    dir_fs = model.FileSystem({SeriesAttributeName.TRUNK: trunk_ds,
                               SeriesAttributeName.BRANCH: branch_ds,
                               SeriesAttributeName.LEAF: leaf_ds})
    return dir_fs


def tau(prefix):
    """ construct the tau filesystem

    layers:
     - trunk (specifiers: [])
     - leaf (specifiers: [tid])

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = dir_.tau_trunk(prefix)
    leaf_ds = dir_.tau_leaf(prefix, root_ds=trunk_ds)

    vma_dfile = file_.vmatrix(FilePrefix.TAU)
    inf_dfile = file_.information(FilePrefix.TAU,
                                  function=info.tau_trunk)
    traj_dfile = file_.trajectory(FilePrefix.TAU)
    trunk_ds.add_data_files({
        FileAttributeName.VMATRIX: vma_dfile,
        FileAttributeName.INFO: inf_dfile,
        FileAttributeName.TRAJ: traj_dfile})

    geom_inf_dfile = file_.information(FilePrefix.GEOM, function=info.run)
    grad_inf_dfile = file_.information(FilePrefix.GRAD, function=info.run)
    hess_inf_dfile = file_.information(FilePrefix.HESS, function=info.run)
    geom_inp_dfile = file_.input_file(FilePrefix.GEOM)
    grad_inp_dfile = file_.input_file(FilePrefix.GRAD)
    hess_inp_dfile = file_.input_file(FilePrefix.HESS)
    ene_dfile = file_.energy(FilePrefix.GEOM)
    geom_dfile = file_.geometry(FilePrefix.GEOM)
    zmat_dfile = file_.zmatrix(FilePrefix.GEOM)
    grad_dfile = file_.gradient(FilePrefix.GRAD)
    hess_dfile = file_.hessian(FilePrefix.HESS)
    hfreq_dfile = file_.harmonic_frequencies(FilePrefix.HESS)
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

    dir_fs = model.FileSystem({SeriesAttributeName.TRUNK: trunk_ds,
                               SeriesAttributeName.LEAF: leaf_ds})
    return dir_fs


def vpt2(prefix):
    """ construct the vpt2 filesystem

    layers:
     - trunk (specifiers: [])

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = dir_.vpt2_trunk(prefix)

    vpt2_inf_dfile = file_.information(FilePrefix.VPT2, function=info.run)
    vpt2_inp_dfile = file_.input_file(FilePrefix.VPT2)
    xmat_dfile = file_.anharmonicity_matrix(FilePrefix.VPT2)

    trunk_ds.add_data_files({
        FileAttributeName.VPT2_INFO: vpt2_inf_dfile,
        FileAttributeName.VPT2_INPUT: vpt2_inp_dfile,
        FileAttributeName.XMAT: xmat_dfile})

    dir_fs = model.FileSystem({SeriesAttributeName.TRUNK: trunk_ds})

    return dir_fs


def reaction(prefix):
    """ construct the reaction filesystem

    layers:
     - trunk (specifiers: [])
     - leaf (specifiers: [rxn_ichs, rxn_chgs, rxn_muls, ts_mul])

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = dir_.reaction_trunk(prefix)
    leaf_ds = dir_.reaction_leaf(prefix, root_ds=trunk_ds)

    dir_fs = model.FileSystem({SeriesAttributeName.TRUNK: trunk_ds,
                               SeriesAttributeName.LEAF: leaf_ds})
    return dir_fs


def ts(prefix):
    """ construct the ts filesystem

    layers:
     - trunk (specifiers: [])

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = dir_.ts_trunk(prefix)

    geom_dfile = file_.geometry(FilePrefix.GEOM)
    ene_dfile = file_.energy(FilePrefix.GEOM)
    zmat_dfile = file_.zmatrix(FilePrefix.GEOM)
    trunk_ds.add_data_files({
        FileAttributeName.ENERGY: ene_dfile,
        FileAttributeName.GEOM: geom_dfile,
        FileAttributeName.ZMAT: zmat_dfile})

    dir_fs = model.FileSystem({SeriesAttributeName.TRUNK: trunk_ds})
    return dir_fs


def direction(prefix):
    """ filesystem object for reaction direction

    layers:
     - leaf (specifiers: [forw])

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    leaf_ds = dir_.direction_leaf(prefix)

    inf_dfile = file_.information(FilePrefix.GEOM, function=info.run)
    inp_dfile = file_.input_file(FilePrefix.GEOM)
    ene_dfile = file_.energy(FilePrefix.GEOM)
    geom_dfile = file_.geometry(FilePrefix.GEOM)
    zmat_dfile = file_.zmatrix(FilePrefix.GEOM)

    leaf_ds.add_data_files({
        FileAttributeName.GEOM_INPUT: inp_dfile,
        FileAttributeName.GEOM_INFO: inf_dfile,
        FileAttributeName.ENERGY: ene_dfile,
        FileAttributeName.GEOM: geom_dfile,
        FileAttributeName.ZMAT: zmat_dfile})

    dir_fs = model.FileSystem({SeriesAttributeName.LEAF: leaf_ds})
    return dir_fs


def run(prefix):
    """ construct the run filesystem [trunk/leaf]

    layers:
     - trunk (specifiers: [])
     - leaf (specifiers: [job])

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = dir_.run_trunk(prefix)
    leaf_ds = dir_.run_leaf(prefix, root_ds=trunk_ds)

    inf_dfile = file_.information(FilePrefix.RUN, function=info.run)
    inp_dfile = file_.input_file(FilePrefix.RUN)
    out_dfile = file_.output_file(FilePrefix.RUN)
    leaf_ds.add_data_files({
        FileAttributeName.INFO: inf_dfile,
        FileAttributeName.INPUT: inp_dfile,
        FileAttributeName.OUTPUT: out_dfile})

    dir_fs = model.FileSystem({SeriesAttributeName.TRUNK: trunk_ds,
                               SeriesAttributeName.LEAF: leaf_ds})
    return dir_fs


def subrun(prefix):
    """ construct the subrun filesystem [leaf]

    layers:
     - leaf (specifiers: [macro_idx, micro_idx])

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    leaf_ds = dir_.subrun_leaf(prefix)

    inf_dfile = file_.information(FilePrefix.RUN, function=info.run)
    inp_dfile = file_.input_file(FilePrefix.RUN)
    out_dfile = file_.output_file(FilePrefix.RUN)
    leaf_ds.add_data_files({
        FileAttributeName.INFO: inf_dfile,
        FileAttributeName.INPUT: inp_dfile,
        FileAttributeName.OUTPUT: out_dfile})

    dir_fs = model.FileSystem({SeriesAttributeName.LEAF: leaf_ds})
    return dir_fs


def build(prefix):
    """ construct the build filesystem

    layers:
     - trunk (specifiers: [head])
     - leaf (specifiers: [head, num])

    :param prefix: sets the path where this filesystem will sit
    :type prefix: str
    """
    trunk_ds = dir_.build_trunk(prefix)
    leaf_ds = dir_.build_leaf(prefix, root_ds=trunk_ds)

    inp_dfile = file_.input_file(FilePrefix.BUILD)
    out_dfile = file_.output_file(FilePrefix.BUILD)
    leaf_ds.add_data_files({
        FileAttributeName.INPUT: inp_dfile,
        FileAttributeName.OUTPUT: out_dfile})

    dir_fs = model.FileSystem({SeriesAttributeName.TRUNK: trunk_ds,
                               SeriesAttributeName.LEAF: leaf_ds})
    return dir_fs


def _process_root_args(root_fs=None, top_ds_name=None):
    if root_fs is not None:
        root_fs = dict(root_fs)
        assert top_ds_name in root_fs
        top_dsdir = root_fs[top_ds_name].dir
    else:
        root_fs = {}
        top_dsdir = None
    return root_fs, top_dsdir
