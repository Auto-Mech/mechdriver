""" utilites """
import autofile

def reaction_energy(save_prefix, rxn_ich, rxn_chg, rxn_mul, method, basis, restrict_open_shell):
    """ reaction energy """
    rct_ichs, prd_ichs = rxn_ich
    rct_chgs, prd_chgs = rxn_chg
    rct_muls, prd_muls = rxn_mul
    rct_enes = reagent_energies(save_prefix, rct_ichs, rct_chgs, rct_muls, method, basis, restrict_open_shell)
    prd_enes = reagent_energies(save_prefix, prd_ichs, prd_chgs, prd_muls, method, basis, restrict_open_shell)
    return sum(prd_enes) - sum(rct_enes)


def reagent_energies(save_prefix, rgt_ichs, rgt_chgs, rgt_muls, method, basis, restrict_open_shell):
    """ reagent energies """
    enes = []
    spc_afs = autofile.fs.species()
    thy_afs = autofile.fs.theory(spc_afs, 'species')
    cnf_afs = autofile.fs.conformer(thy_afs, 'theory')

    for rgt_ich, rgt_chg, rgt_mul in zip(rgt_ichs, rgt_chgs, rgt_muls):
        print(rgt_ich, rgt_mul)
        if restrict_open_shell:
            orb_restr = True
        else:
            orb_restr = (rgt_mul == 1)
        spc_alocs = [rgt_ich, rgt_chg, rgt_mul]
        thy_rlocs = [method, basis, orb_restr]
        thy_alocs = spc_alocs + thy_rlocs
        ene = cnf_afs.conf_trunk.file.energy.read(save_prefix, thy_alocs)
        enes.append(ene)
    return enes

