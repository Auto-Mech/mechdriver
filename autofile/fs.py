""" filesystem library
"""
from autofile.system import model
from autofile.system import series


class AttributeName():
    """ DataFile attribute names """
    REF_TRUNK = 'reference'
    DIR_LEAF = 'direction'
    SPC_TRUNK = 'species_trunk'
    SPC_LEAF = 'species'
    TS_TRUNK = 'ts'
    RXN_TRUNK = 'reaction_trunk'
    RXN_LEAF = 'reaction'
    THY_LEAF = 'theory'
    CNF_TRUNK = 'conf_trunk'
    CNF_LEAF = 'conf'
    SP_TRUNK = 'conf_sp_trunk'
    SP_LEAF = 'conf_sp'
    SCN_TRUNK = 'scan_trunk'
    SCN_BRANCH = 'scan_branch'
    SCN_LEAF = 'scan'
    TAU_TRUNK = 'tau_trunk'
    TAU_LEAF = 'tau'
    SPT_TRUNK = 'tau_sp_trunk'
    SPT_LEAF = 'tau_sp'
    RUN_TRUNK = 'run_trunk'
    RUN_LEAF = 'run'
    BUILD_TRUNK = 'build_trunk'
    BUILD_LEAF = 'build'


def empty():
    """ create an empty filesystem
    """
    return model.FileSystem({})


def reference(root_fs=None, top_ds_name=None, name_prefix=''):
    """ construct the reference filesystem
    """
    root_fs, top_dsdir = _process_root_args(root_fs, top_ds_name)

    ref_trunk_ds = series.reference_trunk(root_dsdir=top_dsdir)

    ref_fs = model.FileSystem({
        (name_prefix + AttributeName.REF_TRUNK): ref_trunk_ds,
    })
    ref_fs.update(root_fs)
    return ref_fs


def direction(root_fs=None, top_ds_name=None, name_prefix=''):
    """ construct the direction filesystem
    """
    root_fs, top_dsdir = _process_root_args(root_fs, top_ds_name)

    dir_trunk_ds = series.direction_leaf(root_dsdir=top_dsdir)

    dir_fs = model.FileSystem({
        (name_prefix + AttributeName.DIR_LEAF): dir_trunk_ds,
    })
    dir_fs.update(root_fs)
    return dir_fs


def species(root_fs=None, top_ds_name=None, name_prefix=''):
    """ construct the species filesystem
    """
    root_fs, top_dsdir = _process_root_args(root_fs, top_ds_name)

    spc_trunk_ds = series.species_trunk(root_dsdir=top_dsdir)
    spc_leaf_ds = series.species_leaf(root_dsdir=spc_trunk_ds.dir)

    spc_fs = model.FileSystem({
        (name_prefix + AttributeName.SPC_TRUNK): spc_trunk_ds,
        (name_prefix + AttributeName.SPC_LEAF): spc_leaf_ds,
    })
    spc_fs.update(root_fs)
    return spc_fs


def ts(root_fs=None, top_ds_name=None, name_prefix=''):
    """ construct the species filesystem
    """
    root_fs, top_dsdir = _process_root_args(root_fs, top_ds_name)

    ts_trunk_ds = series.ts_trunk(root_dsdir=top_dsdir)

    ts_fs = model.FileSystem({
        (name_prefix + AttributeName.TS_TRUNK): ts_trunk_ds,
    })
    ts_fs.update(root_fs)
    return ts_fs


def reaction(root_fs=None, top_ds_name=None, name_prefix=''):
    """ construct the reaction filesystem
    """
    root_fs, top_dsdir = _process_root_args(root_fs, top_ds_name)

    rxn_trunk_ds = series.reaction_trunk(root_dsdir=top_dsdir)
    rxn_leaf_ds = series.reaction_leaf(root_dsdir=rxn_trunk_ds.dir)

    rxn_fs = model.FileSystem({
        (name_prefix + AttributeName.RXN_TRUNK): rxn_trunk_ds,
        (name_prefix + AttributeName.RXN_LEAF): rxn_leaf_ds,
    })
    rxn_fs.update(root_fs)
    return rxn_fs


def theory(root_fs=None, top_ds_name=None, name_prefix=''):
    """ construct the theory filesystem
    """
    root_fs, top_dsdir = _process_root_args(root_fs, top_ds_name)

    thy_leaf_ds = series.theory_leaf(root_dsdir=top_dsdir)

    thy_fs = model.FileSystem({
        (name_prefix + AttributeName.THY_LEAF): thy_leaf_ds,
    })
    thy_fs.update(root_fs)
    return thy_fs


def conformer(root_fs=None, top_ds_name=None, name_prefix=''):
    """ construct the conformer filesystem
    """
    root_fs, top_dsdir = _process_root_args(root_fs, top_ds_name)

    cnf_trunk_ds = series.conformer_trunk(root_dsdir=top_dsdir)
    cnf_leaf_ds = series.conformer_leaf(root_dsdir=cnf_trunk_ds.dir)

    cnf_fs = model.FileSystem({
        (name_prefix + AttributeName.CNF_TRUNK): cnf_trunk_ds,
        (name_prefix + AttributeName.CNF_LEAF): cnf_leaf_ds,
    })
    cnf_fs.update(root_fs)
    return cnf_fs


def tau(root_fs=None, top_ds_name=None, name_prefix=''):
    """ construct the tau filesystem
    """
    root_fs, top_dsdir = _process_root_args(root_fs, top_ds_name)

    tau_trunk_ds = series.tau_trunk(root_dsdir=top_dsdir)
    tau_leaf_ds = series.tau_leaf(root_dsdir=tau_trunk_ds.dir)

    tau_fs = model.FileSystem({
        (name_prefix + AttributeName.TAU_TRUNK): tau_trunk_ds,
        (name_prefix + AttributeName.TAU_LEAF): tau_leaf_ds,
    })
    tau_fs.update(root_fs)
    return tau_fs


def single_point(root_fs=None, top_ds_name=None, name_prefix=''):
    """ construct the single-point filesystem
    """
    root_fs, top_dsdir = _process_root_args(root_fs, top_ds_name)

    sp_trunk_ds = series.single_point_trunk(root_dsdir=top_dsdir)
    sp_leaf_ds = series.single_point_leaf(root_dsdir=sp_trunk_ds.dir)

    sp_fs = model.FileSystem({
        (name_prefix + AttributeName.SP_TRUNK): sp_trunk_ds,
        (name_prefix + AttributeName.SP_LEAF): sp_leaf_ds,
    })
    sp_fs.update(root_fs)
    return sp_fs


def scan(root_fs=None, top_ds_name=None, name_prefix=''):
    """ construct the scan filesystem
    """
    root_fs, top_dsdir = _process_root_args(root_fs, top_ds_name)

    scn_trunk_ds = series.scan_trunk(root_dsdir=top_dsdir)
    scn_branch_ds = series.scan_branch(root_dsdir=scn_trunk_ds.dir)
    scn_leaf_ds = series.scan_leaf(root_dsdir=scn_branch_ds.dir)

    scn_fs = model.FileSystem({
        (name_prefix + AttributeName.SCN_TRUNK): scn_trunk_ds,
        (name_prefix + AttributeName.SCN_BRANCH): scn_branch_ds,
        (name_prefix + AttributeName.SCN_LEAF): scn_leaf_ds,
    })
    scn_fs.update(root_fs)
    return scn_fs


def run(root_fs=None, top_ds_name=None, name_prefix=''):
    """ construct the run filesystem
    """
    root_fs, top_dsdir = _process_root_args(root_fs, top_ds_name)

    run_trunk_ds = series.run_trunk(root_dsdir=top_dsdir)
    run_leaf_ds = series.run_leaf(root_dsdir=run_trunk_ds.dir)

    run_fs = model.FileSystem({
        (name_prefix + AttributeName.RUN_TRUNK): run_trunk_ds,
        (name_prefix + AttributeName.RUN_LEAF): run_leaf_ds,
    })
    run_fs.update(root_fs)
    return run_fs


def build(root_fs=None, top_ds_name=None, name_prefix=''):
    """ construct the run filesystem
    """
    root_fs, top_dsdir = _process_root_args(root_fs, top_ds_name)

    build_trunk_ds = series.build_trunk(root_dsdir=top_dsdir)
    build_leaf_ds = series.build_leaf(root_dsdir=build_trunk_ds.dir)

    build_fs = model.FileSystem({
        (name_prefix + AttributeName.BUILD_TRUNK): build_trunk_ds,
        (name_prefix + AttributeName.BUILD_LEAF): build_leaf_ds,
    })
    build_fs.update(root_fs)
    return build_fs


def _process_root_args(root_fs=None, top_ds_name=None):
    if root_fs is not None:
        root_fs = dict(root_fs)
        assert top_ds_name in root_fs
        top_dsdir = root_fs[top_ds_name].dir
    else:
        root_fs = {}
        top_dsdir = None
    return root_fs, top_dsdir
