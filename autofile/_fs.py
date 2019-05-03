""" creates the filesystem
"""
from autofile.system import model
from autofile.system import series


class AttributeName():
    """ DataFile attribute names """
    SPC_TRUNK = 'species_trunk'
    SPC_LEAF = 'species'
    THY_LEAF = 'theory'
    CNF_TRUNK = 'conformer_trunk'
    CNF_LEAF = 'conformer'
    CNF_RUN_LEAF = 'conformer_run'
    SCN_TRUNK = 'scan_trunk'
    SCN_BRANCH = 'scan_branch'
    SCN_LEAF = 'scan'
    SCN_RUN_LEAF = 'scan_run'


SPC_TRUNK_DS = series.species_trunk()
SPC_LEAF_DS = series.species_leaf(root_dsdir=SPC_TRUNK_DS.dir)
THY_LEAF_DS = series.theory_leaf(root_dsdir=SPC_LEAF_DS.dir)
CNF_TRUNK_DS = series.conformer_trunk(root_dsdir=THY_LEAF_DS.dir)
CNF_LEAF_DS = series.conformer_leaf(root_dsdir=CNF_TRUNK_DS.dir)
CNF_RUN_TRUNK_DS = series.run_trunk(root_dsdir=CNF_LEAF_DS.dir)
CNF_RUN_LEAF_DS = series.run_leaf(root_dsdir=CNF_RUN_TRUNK_DS.dir)
SCN_TRUNK_DS = series.scan_trunk(root_dsdir=CNF_LEAF_DS.dir)
SCN_BRANCH_DS = series.scan_branch(root_dsdir=SCN_TRUNK_DS.dir)
SCN_LEAF_DS = series.scan_leaf(root_dsdir=SCN_BRANCH_DS.dir)
SCN_RUN_TRUNK_DS = series.run_trunk(root_dsdir=SCN_LEAF_DS.dir)
SCN_RUN_LEAF_DS = series.run_leaf(root_dsdir=SCN_RUN_TRUNK_DS.dir)

FS_ = model.FileSystem({
    AttributeName.SPC_TRUNK: SPC_TRUNK_DS,
    AttributeName.SPC_LEAF: SPC_LEAF_DS,
    AttributeName.THY_LEAF: THY_LEAF_DS,
    AttributeName.CNF_TRUNK: CNF_TRUNK_DS,
    AttributeName.CNF_LEAF: CNF_LEAF_DS,
    AttributeName.CNF_RUN_LEAF: CNF_RUN_LEAF_DS,
    AttributeName.SCN_TRUNK: SCN_TRUNK_DS,
    AttributeName.SCN_BRANCH: SCN_BRANCH_DS,
    AttributeName.SCN_LEAF: SCN_LEAF_DS,
    AttributeName.SCN_RUN_LEAF: SCN_RUN_LEAF_DS,
})
