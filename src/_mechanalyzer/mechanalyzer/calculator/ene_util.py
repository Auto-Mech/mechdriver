""" functions for energy manipulation
"""

from phydat import phycon

def max_en_auto(n_atoms, ene_bw, ref_ene=0, T=2500):
    """ Determines automatically the max energy of microcanonical output and max
        energy stored in products for appropriate ped and hoten output writing
        :param n_atoms: number of atoms involved in the bimol reaction
        :type n_atoms: int
        :param ene_bw: backward energy barrier of the reaction
                       > 0 for exothermic reactions
                       < 0 for enothermic reactions
        :type ene_bw: float
        :param ref_ene: reference energy (bimol prods energy in ped)
        :type ref_ene: float
        :param T: reference temperature
        :type T: float
        :return max_ene: maximum energy written in mess ped/hoten output

    """
    # determine average boltzmann energy for the TS at 2500 K
    # NB equipartition theorem - no trasl bc preserved in rxn
    # assume non-linear
    boltz_ene_T = phycon.RC_KCAL*T*(3*n_atoms-7+3/2)
    sigma = 0.87+0.04*(boltz_ene_T + ene_bw)  # from Danilack 2020
    print('debug ene_util line 30 energies: ref {}, boltz {}, barrier {}, sigma {} \n'.format(
        ref_ene, boltz_ene_T, ene_bw, sigma))
    max_ene = ref_ene + boltz_ene_T + ene_bw + 4*sigma

    return max_ene

