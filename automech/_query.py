"""Functions for querying the filesystem."""

import itertools

import autofile
import automol


# Paths from root
def reaction_path_from_root(root_path: str, rsmis, psmis) -> str:
    """Get reaction file system path from the root path, given SMILES

    :param root_path: The root path of the file system
    :param rsmis: Reactant SMILES
    :param psmis: Product SMILES
    :returns: The file system path
    """
    rchis = list(map(automol.smiles.chi, rsmis))
    pchis = list(map(automol.smiles.chi, psmis))
    rxn_fs = autofile.fs.reaction(root_path)
    for locs in rxn_fs[-1].existing():
        (rchis0, pchis0), *_ = locs
        for rchis1, pchis1 in itertools.permutations([rchis0, pchis0]):
            if rchis1 == rchis and pchis1 == pchis:
                return rxn_fs[-1].path(locs)


def transition_state_paths_from_root(
    root_path: str, rsmis, psmis, method, basis
) -> list[str]:
    """Get transition state paths from the root path, given SMILES, method, and basis

    :param root_path: The root path of the file system
    :param rsmis: Reactant SMILES
    :param psmis: Product SMILES
    :param method: The electronic structure method
    :param basis: The electronic structure basis set
    :return: The transition state paths
    """
    rxn_path = reaction_path_from_root(root_path, rsmis=rsmis, psmis=psmis)
    thy_path = theory_path_from_prefix(rxn_path, method, basis)
    ts_fs = autofile.fs.transition_state(thy_path)
    return [ts_fs[-1].path(locs) for locs in ts_fs[-1].existing()]


# Paths from prefix
def theory_path_from_prefix(prefix: str, method: str, basis: str) -> str:
    """Get theory file system path from its prefix, given method and basis

    :param prefix: The prefix of the theory file system
    :param method: The electronic structure method
    :param basis: The electronic structure basis set
    :return: The file system path
    """
    thy_fs = autofile.fs.theory(prefix)
    method = method.lower()
    basis = basis.lower()
    for locs in thy_fs[-1].existing():
        method0, basis0, *_ = locs
        if method == method0 and basis == basis0:
            return thy_fs[-1].path(locs)


def conformers_locators_from_prefix(prefix: str) -> list[object]:
    """Get locators for a set of conformers from their prefix

    :param prefix: The prefix of the conformer file system
    :return: The locator values
    """
    cnf_fs = autofile.fs.conformer(prefix)
    return list(cnf_fs[-1].existing())


def visualize_conformers_from_prefix(prefix: str):
    """Visualize the conformers at a given prefix

    :param prefix: The prefix of the conformer file system
    """
    cnf_fs = autofile.fs.conformer(prefix)
    for cnf_locs in conformers_locators_from_prefix(prefix):
        geo = cnf_fs[-1].file.geometry.read(cnf_locs)
        hess = cnf_fs[-1].file.hessian.read(cnf_locs)
        automol.geom.display(geo)
        print(hess)
