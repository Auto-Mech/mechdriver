""" utilites """
import autofile

def reaction_energy(save_prefix, rxn_ich, rxn_chg, rxn_mul, meth, bas, restrict_open_shell):
    """ reaction energy """
    rct_ichs_, prd_ichs_ = rxn_ich
    rct_chgs_, prd_chgs_ = rxn_chg
    rct_muls_, prd_muls_ = rxn_mul
    rct_enes = reagent_energies(save_prefix, rct_ichs_, rct_chgs_, rct_muls_, meth, bas, restrict_open_shell)
    prd_enes = reagent_energies(save_prefix, prd_ichs_, prd_chgs_, prd_muls_, meth, bas, restrict_open_shell)
    return sum(prd_enes) - sum(rct_enes)


def reagent_energies(save_prefix, rgt_ichs, rgt_chgs, rgt_muls, meth, bas, restrict_open_shell):
    """ reagent energies """
    enes_ = []
    spc_afs = autofile.fs.species()
    thy_afs = autofile.fs.theory(spc_afs, 'species')
    cnf_afs = autofile.fs.conformer(thy_afs, 'theory')

    for rgt_ich, rgt_chg, rgt_mul in zip(rgt_ichs, rgt_chgs, rgt_muls):
        print(rgt_ich, rgt_mul)
        if restrict_open_shell:
            orb_restr_ = True
        else:
            orb_restr_ = (rgt_mul == 1)
        spc_alocs_ = [rgt_ich, rgt_chg, rgt_mul]
        thy_rlocs_ = [meth, bas, orb_restr_]
        thy_alocs_ = spc_alocs_ + thy_rlocs_
        ene_ = cnf_afs.conf_trunk.file.energy.read(save_prefix, thy_alocs_)
        enes_.append(ene_)
    return enes_

