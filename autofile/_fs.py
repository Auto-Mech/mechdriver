""" creates the filesystem
"""
from autofile.system import model
from autofile.system import series


class AttributeName():
    """ DataFile attribute names """
    SPC_TRUNK = 'species_trunk'
    SPC_LEAF = 'species'
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


SPC_TRUNK_DS = series.species_trunk()
SPC_LEAF_DS = series.species_leaf(root_dsdir=SPC_TRUNK_DS.dir)
THY_LEAF_DS = series.theory_leaf(root_dsdir=SPC_LEAF_DS.dir)
CNF_TRUNK_DS = series.conformer_trunk(root_dsdir=THY_LEAF_DS.dir)
CNF_LEAF_DS = series.conformer_leaf(root_dsdir=CNF_TRUNK_DS.dir)
CNF_RUN_TRUNK_DS = series.run_trunk(root_dsdir=CNF_LEAF_DS.dir)
CNF_RUN_LEAF_DS = series.run_leaf(root_dsdir=CNF_RUN_TRUNK_DS.dir)
TAU_TRUNK_DS = series.tau_trunk(root_dsdir=THY_LEAF_DS.dir)
TAU_LEAF_DS = series.tau_leaf(root_dsdir=TAU_TRUNK_DS.dir)
TAU_RUN_TRUNK_DS = series.run_trunk(root_dsdir=TAU_LEAF_DS.dir)
TAU_RUN_LEAF_DS = series.run_leaf(root_dsdir=TAU_RUN_TRUNK_DS.dir)

# subdirectories of conformer
# 1. single-point
SP_TRUNK_DS = series.single_point_trunk(root_dsdir=CNF_LEAF_DS.dir)
SP_LEAF_DS = series.single_point_leaf(root_dsdir=SP_TRUNK_DS.dir)
# 2. scan
SCN_TRUNK_DS = series.scan_trunk(root_dsdir=CNF_LEAF_DS.dir)
SCN_BRANCH_DS = series.scan_branch(root_dsdir=SCN_TRUNK_DS.dir)
SCN_LEAF_DS = series.scan_leaf(root_dsdir=SCN_BRANCH_DS.dir)
SCN_RUN_TRUNK_DS = series.run_trunk(root_dsdir=SCN_LEAF_DS.dir)
SCN_RUN_LEAF_DS = series.run_leaf(root_dsdir=SCN_RUN_TRUNK_DS.dir)
# subdirectories of tau
# 1. single-point
SPT_TRUNK_DS = series.single_point_trunk(root_dsdir=TAU_LEAF_DS.dir)
SPT_LEAF_DS = series.single_point_leaf(root_dsdir=SPT_TRUNK_DS.dir)

FS_ = model.FileSystem({
    AttributeName.SPC_TRUNK: SPC_TRUNK_DS,
    AttributeName.SPC_LEAF: SPC_LEAF_DS,
    AttributeName.THY_LEAF: THY_LEAF_DS,
    AttributeName.CNF_TRUNK: CNF_TRUNK_DS,
    AttributeName.CNF_LEAF: CNF_LEAF_DS,
    AttributeName.CNF_RUN_LEAF: CNF_RUN_LEAF_DS,
    AttributeName.SP_TRUNK: SP_TRUNK_DS,
    AttributeName.SP_LEAF: SP_LEAF_DS,
    AttributeName.SCN_TRUNK: SCN_TRUNK_DS,
    AttributeName.SCN_BRANCH: SCN_BRANCH_DS,
    AttributeName.SCN_LEAF: SCN_LEAF_DS,
    AttributeName.SCN_RUN_LEAF: SCN_RUN_LEAF_DS,
    AttributeName.TAU_TRUNK: TAU_TRUNK_DS,
    AttributeName.TAU_LEAF: TAU_LEAF_DS,
    AttributeName.TAU_RUN_LEAF: TAU_RUN_LEAF_DS,
    AttributeName.SPT_TRUNK: SPT_TRUNK_DS,
    AttributeName.SPT_LEAF: SPT_LEAF_DS,
})
