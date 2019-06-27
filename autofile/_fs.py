""" creates the species filesystem
"""
from autofile.system import model
from autofile.system import series


class AttributeName():
    """ DataFile attribute names """
    SPC_TRUNK = 'species_trunk'
    SPC_LEAF = 'species'
    RXN_TRUNK = 'reaction_trunk'
    RXN_LEAF = 'reaction'
    THY_LEAF = 'theory'
    CNF_TRUNK = 'conf_trunk'
    CNF_LEAF = 'conf'
    CNF_RUN_LEAF = 'conf_run'
    SP_TRUNK = 'conf_sp_trunk'
    SP_LEAF = 'conf_sp'
    SCN_TRUNK = 'scan_trunk'
    SCN_BRANCH = 'scan_branch'
    SCN_LEAF = 'scan'
    SCN_RUN_LEAF = 'scan_run'
    TAU_TRUNK = 'tau_trunk'
    TAU_LEAF = 'tau'
    TAU_RUN_LEAF = 'tau_run'
    SPT_TRUNK = 'tau_sp_trunk'
    SPT_LEAF = 'tau_sp'


def species_filesystem():
    """ construct the species filesystem
    """
    spc_trunk_ds = series.species_trunk()
    spc_leaf_ds = series.species_leaf(root_dsdir=spc_trunk_ds.dir)
    thy_leaf_ds = series.theory_leaf(root_dsdir=spc_leaf_ds.dir)
    cnf_trunk_ds = series.conformer_trunk(root_dsdir=thy_leaf_ds.dir)
    cnf_leaf_ds = series.conformer_leaf(root_dsdir=cnf_trunk_ds.dir)
    cnf_run_trunk_ds = series.run_trunk(root_dsdir=cnf_leaf_ds.dir)
    cnf_run_leaf_ds = series.run_leaf(root_dsdir=cnf_run_trunk_ds.dir)
    tau_trunk_ds = series.tau_trunk(root_dsdir=thy_leaf_ds.dir)
    tau_leaf_ds = series.tau_leaf(root_dsdir=tau_trunk_ds.dir)
    tau_run_trunk_ds = series.run_trunk(root_dsdir=tau_leaf_ds.dir)
    tau_run_leaf_ds = series.run_leaf(root_dsdir=tau_run_trunk_ds.dir)

    # subdirectories of conformer
    # 1. single-point
    sp_trunk_ds = series.single_point_trunk(root_dsdir=cnf_leaf_ds.dir)
    sp_leaf_ds = series.single_point_leaf(root_dsdir=sp_trunk_ds.dir)
    # 2. scan
    scn_trunk_ds = series.scan_trunk(root_dsdir=cnf_leaf_ds.dir)
    scn_branch_ds = series.scan_branch(root_dsdir=scn_trunk_ds.dir)
    scn_leaf_ds = series.scan_leaf(root_dsdir=scn_branch_ds.dir)
    scn_run_trunk_ds = series.run_trunk(root_dsdir=scn_leaf_ds.dir)
    scn_run_leaf_ds = series.run_leaf(root_dsdir=scn_run_trunk_ds.dir)
    # subdirectories of tau
    # 1. single-point
    spt_trunk_ds = series.single_point_trunk(root_dsdir=tau_leaf_ds.dir)
    spt_leaf_ds = series.single_point_leaf(root_dsdir=spt_trunk_ds.dir)

    spc_fs = model.FileSystem({
        AttributeName.SPC_TRUNK: spc_trunk_ds,
        AttributeName.SPC_LEAF: spc_leaf_ds,
        AttributeName.THY_LEAF: thy_leaf_ds,
        AttributeName.CNF_TRUNK: cnf_trunk_ds,
        AttributeName.CNF_LEAF: cnf_leaf_ds,
        AttributeName.CNF_RUN_LEAF: cnf_run_leaf_ds,
        AttributeName.SP_TRUNK: sp_trunk_ds,
        AttributeName.SP_LEAF: sp_leaf_ds,
        AttributeName.SCN_TRUNK: scn_trunk_ds,
        AttributeName.SCN_BRANCH: scn_branch_ds,
        AttributeName.SCN_LEAF: scn_leaf_ds,
        AttributeName.SCN_RUN_LEAF: scn_run_leaf_ds,
        AttributeName.TAU_TRUNK: tau_trunk_ds,
        AttributeName.TAU_LEAF: tau_leaf_ds,
        AttributeName.TAU_RUN_LEAF: tau_run_leaf_ds,
        AttributeName.SPT_TRUNK: spt_trunk_ds,
        AttributeName.SPT_LEAF: spt_leaf_ds,
    })

    return spc_fs


def reaction_filesystem():
    """ construct the reaction filesystem
    """
    rxn_trunk_ds = series.reaction_trunk()
    rxn_leaf_ds = series.reaction_leaf(root_dsdir=rxn_trunk_ds.dir)
    thy_leaf_ds = series.theory_leaf(root_dsdir=rxn_leaf_ds.dir)
    cnf_trunk_ds = series.conformer_trunk(root_dsdir=thy_leaf_ds.dir)
    cnf_leaf_ds = series.conformer_leaf(root_dsdir=cnf_trunk_ds.dir)
    cnf_run_trunk_ds = series.run_trunk(root_dsdir=cnf_leaf_ds.dir)
    cnf_run_leaf_ds = series.run_leaf(root_dsdir=cnf_run_trunk_ds.dir)

    # (gridopt) scan directories
    scn_trunk_ds = series.scan_trunk(root_dsdir=cnf_leaf_ds.dir)
    scn_branch_ds = series.scan_branch(root_dsdir=scn_trunk_ds.dir)
    scn_leaf_ds = series.scan_leaf(root_dsdir=scn_branch_ds.dir)
    scn_run_trunk_ds = series.run_trunk(root_dsdir=scn_leaf_ds.dir)
    scn_run_leaf_ds = series.run_leaf(root_dsdir=scn_run_trunk_ds.dir)

    rxn_fs = model.FileSystem({
        AttributeName.RXN_TRUNK: rxn_trunk_ds,
        AttributeName.RXN_LEAF: rxn_leaf_ds,
        AttributeName.THY_LEAF: thy_leaf_ds,
        AttributeName.CNF_TRUNK: cnf_trunk_ds,
        AttributeName.CNF_LEAF: cnf_leaf_ds,
        AttributeName.CNF_RUN_LEAF: cnf_run_leaf_ds,
        AttributeName.SCN_TRUNK: scn_trunk_ds,
        AttributeName.SCN_BRANCH: scn_branch_ds,
        AttributeName.SCN_LEAF: scn_leaf_ds,
        AttributeName.SCN_RUN_LEAF: scn_run_leaf_ds,
    })

    return rxn_fs
