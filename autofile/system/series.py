""" DataSeries
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


class FileAttributeName():
    """ DataFile attribute names """
    INFO = 'info'
    INPUT = 'input'
    OUTPUT = 'output'
    VMATRIX = 'vmatrix'
    GEOM_INFO = 'geometry_info'
    GRAD_INFO = 'gradient_info'
    HESS_INFO = 'hessian_info'
    GEOM_INPUT = 'geometry_input'
    GRAD_INPUT = 'gradient_input'
    HESS_INPUT = 'hessian_input'
    ENERGY = 'energy'
    GEOM = 'geometry'
    ZMAT = 'zmatrix'
    GRAD = 'gradient'
    HESS = 'hessian'
    HFREQ = 'harmonic_frequencies'
    TRAJ = 'trajectory'


def species_trunk(root_dsdir=None):
    """ species trunk DataSeries
    """
    dsdir = dir_.species_trunk(root_dsdir)
    return model.DataSeries(dsdir=dsdir)


def ts_trunk(root_dsdir=None):
    """ ts trunk DataSeries
    """
    dsdir = dir_.ts_trunk(root_dsdir)
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

    dseries = model.DataSeries(
        dsdir=dsdir,
        dfile_dct={
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
    return dseries


def direction_leaf(root_dsdir=None):
    """ direction leaf DataSeries
    """
    dsdir = dir_.direction_leaf(root_dsdir)
    inf_dfile = file_.information(FilePrefix.GEOM, function=info.run)
    inp_dfile = file_.input_file(FilePrefix.GEOM)
    ene_dfile = file_.energy(FilePrefix.GEOM)
    geom_dfile = file_.geometry(FilePrefix.GEOM)
    zmat_dfile = file_.zmatrix(FilePrefix.GEOM)

    dseries = model.DataSeries(
        dsdir=dsdir,
        dfile_dct={
            FileAttributeName.GEOM_INPUT: inp_dfile,
            FileAttributeName.GEOM_INFO: inf_dfile,
            FileAttributeName.ENERGY: ene_dfile,
            FileAttributeName.GEOM: geom_dfile,
            FileAttributeName.ZMAT: zmat_dfile})
    return dseries


def species_leaf(root_dsdir=None):
    """ species leaf DataSeries
    """
    dsdir = dir_.species_leaf(root_dsdir)
    return model.DataSeries(dsdir=dsdir)


def reaction_trunk(root_dsdir=None):
    """ reaction trunk DataSeries
    """
    dsdir = dir_.reaction_trunk(root_dsdir)
    return model.DataSeries(dsdir=dsdir)


def reaction_leaf(root_dsdir=None):
    """ reaction leaf DataSeries
    """
    dsdir = dir_.reaction_leaf(root_dsdir)
    return model.DataSeries(dsdir=dsdir)


def theory_leaf(root_dsdir=None):
    """ theory leaf DataSeries
    """
    dsdir = dir_.theory_leaf(root_dsdir)
    geom_dfile = file_.geometry(FilePrefix.GEOM)
    ene_dfile = file_.energy(FilePrefix.GEOM)
    zmat_dfile = file_.zmatrix(FilePrefix.GEOM)

    dseries = model.DataSeries(
        dsdir=dsdir,
        dfile_dct={
            FileAttributeName.ENERGY: ene_dfile,
            FileAttributeName.GEOM: geom_dfile,
            FileAttributeName.ZMAT: zmat_dfile})
    return dseries


def run_trunk(root_dsdir=None):
    """ run trunk DataSeries
    """
    dsdir = dir_.run_trunk(root_dsdir)
    return model.DataSeries(dsdir=dsdir)


def run_leaf(root_dsdir=None):
    """ run leaf DataSeries
    """
    dsdir = dir_.run_leaf(root_dsdir)
    inf_dfile = file_.information(FilePrefix.RUN, function=info.run)
    inp_dfile = file_.input_file(FilePrefix.RUN)
    out_dfile = file_.output_file(FilePrefix.RUN)
    dseries = model.DataSeries(
        dsdir=dsdir,
        dfile_dct={
            FileAttributeName.INFO: inf_dfile,
            FileAttributeName.INPUT: inp_dfile,
            FileAttributeName.OUTPUT: out_dfile})
    return dseries


def subrun_leaf(root_dsdir=None):
    """ subrun leaf DataSeries
    """
    dsdir = dir_.subrun_leaf(root_dsdir)
    inf_dfile = file_.information(FilePrefix.RUN, function=info.run)
    inp_dfile = file_.input_file(FilePrefix.RUN)
    out_dfile = file_.output_file(FilePrefix.RUN)
    dseries = model.DataSeries(
        dsdir=dsdir,
        dfile_dct={
            FileAttributeName.INFO: inf_dfile,
            FileAttributeName.INPUT: inp_dfile,
            FileAttributeName.OUTPUT: out_dfile})
    return dseries


def conformer_trunk(root_dsdir=None):
    """ conformer trunk DataSeries
    """
    dsdir = dir_.conformer_trunk(root_dsdir)
    min_ene_dfile = file_.energy(FilePrefix.MIN)
    vma_dfile = file_.vmatrix(FilePrefix.CONF)
    inf_dfile = file_.information(FilePrefix.CONF,
                                  function=info.conformer_trunk)
    traj_dfile = file_.trajectory(FilePrefix.CONF)

    dseries = model.DataSeries(
        dsdir=dsdir,
        dfile_dct={
            FileAttributeName.VMATRIX: vma_dfile,
            FileAttributeName.INFO: inf_dfile,
            FileAttributeName.ENERGY: min_ene_dfile,
            FileAttributeName.TRAJ: traj_dfile})
    return dseries


def conformer_leaf(root_dsdir=None):
    """ conformer leaf DataSeries
    """
    dsdir = dir_.conformer_leaf(root_dsdir)
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

    dseries = model.DataSeries(
        dsdir=dsdir,
        dfile_dct={
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
    return dseries


def single_point_trunk(root_dsdir=None):
    """ single point trunk DataSeries
    """
    dsdir = dir_.single_point_trunk(root_dsdir)
    return model.DataSeries(dsdir=dsdir)


def single_point_leaf(root_dsdir=None):
    """ single_point leaf DataSeries
    """
    dsdir = dir_.single_point_leaf(root_dsdir)
    ene_dfile = file_.energy(FilePrefix.SP)
    dseries = model.DataSeries(
        dsdir=dsdir,
        dfile_dct={
            FileAttributeName.ENERGY: ene_dfile})
    return dseries


def scan_trunk(root_dsdir=None):
    """ scan trunk DataSeries
    """
    dsdir = dir_.scan_trunk(root_dsdir)
    vma_dfile = file_.vmatrix(FilePrefix.SCAN)

    dseries = model.DataSeries(
        dsdir=dsdir,
        dfile_dct={
            FileAttributeName.VMATRIX: vma_dfile})
    return dseries


def scan_branch(root_dsdir=None):
    """ scan branch DataSeries
    """
    dsdir = dir_.scan_branch(root_dsdir)
    inf_dfile = file_.information(FilePrefix.SCAN, function=info.scan_branch)
    traj_dfile = file_.trajectory(FilePrefix.SCAN)

    dseries = model.DataSeries(
        dsdir=dsdir,
        dfile_dct={
            FileAttributeName.INFO: inf_dfile,
            FileAttributeName.TRAJ: traj_dfile})
    return dseries


def scan_leaf(root_dsdir=None):
    """ scan leaf DataSeries
    """
    dsdir = dir_.scan_leaf(root_dsdir)
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

    dseries = model.DataSeries(
        dsdir=dsdir,
        dfile_dct={
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
    return dseries


def tau_trunk(root_dsdir=None):
    """ tau trunk DataSeries
    """
    dsdir = dir_.tau_trunk(root_dsdir)
    vma_dfile = file_.vmatrix(FilePrefix.TAU)
    inf_dfile = file_.information(FilePrefix.TAU,
                                  function=info.tau_trunk)
    traj_dfile = file_.trajectory(FilePrefix.TAU)

    dseries = model.DataSeries(
        dsdir=dsdir,
        dfile_dct={
            FileAttributeName.VMATRIX: vma_dfile,
            FileAttributeName.INFO: inf_dfile,
            FileAttributeName.TRAJ: traj_dfile})
    return dseries


def tau_leaf(root_dsdir=None):
    """ tau leaf DataSeries
    """
    dsdir = dir_.tau_leaf(root_dsdir)
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

    dseries = model.DataSeries(
        dsdir=dsdir,
        dfile_dct={
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
    return dseries


def build_trunk(root_dsdir=None):
    """ build trunk DataSeries
    """
    dsdir = dir_.build_trunk(root_dsdir)
    return model.DataSeries(dsdir=dsdir)


def build_leaf(root_dsdir=None):
    """ build leaf DataSeries
    """
    dsdir = dir_.build_leaf(root_dsdir)
    inp_dfile = file_.input_file(FilePrefix.BUILD)
    out_dfile = file_.output_file(FilePrefix.BUILD)
    dseries = model.DataSeries(
        dsdir=dsdir,
        dfile_dct={
            FileAttributeName.INPUT: inp_dfile,
            FileAttributeName.OUTPUT: out_dfile})
    return dseries
