""" calculates energies for channels
"""

from phydat import phycon
from mechlib.amech_io import printer as ioprinter
from mechroutines.pf import thermo as thmroutines


# Functions to handle energies for a channel
def set_reference_ene(rxn_lst, spc_dct,
                      pes_model_dct_i, spc_model_dct_i,
                      run_prefix, save_prefix, ref_idx=0):
    """ Sets the reference species for the PES for which all energies
        are scaled relative to.
    """

    # Set the index for the reference species, right now defualt to 1st spc
    ref_rxn = rxn_lst[ref_idx]

    _, (ref_rgts, _) = ref_rxn

    ioprinter.info_message(
        'Determining the reference energy for PES...', newline=1)
    ioprinter.info_message(
        ' - Reference species assumed to be the',
        ' first set of reactants on PES: {}'.format('+'.join(ref_rgts)))

    # Get the model for the first reference species
    ref_scheme = pes_model_dct_i['therm_fit']['ref_scheme']
    ref_enes = pes_model_dct_i['therm_fit']['ref_enes']

    ref_ene_level = spc_model_dct_i['ene']['lvl1'][0]
    ioprinter.info_message(
        ' - Energy Level for Reference Species: {}'.format(ref_ene_level))

    # Get the elec+zpe energy for the reference species
    ioprinter.info_message('')
    hf0k = 0.0
    for rgt in ref_rgts:

        ioprinter.info_message(' - Calculating energy for {}...'.format(rgt))
        basis_dct, uniref_dct = thmroutines.basis.prepare_refs(
            ref_scheme, spc_dct, [[rgt, None]], run_prefix, save_prefix)
        spc_basis, coeff_basis = basis_dct[rgt]

        # Build filesystem
        ene_spc, ene_basis = thmroutines.basis.basis_energy(
            rgt, spc_basis, uniref_dct, spc_dct,
            spc_model_dct_i, run_prefix, save_prefix)

        # Calcualte the total energy
        hf0k += thmroutines.heatform.calc_hform_0k(
            ene_spc, ene_basis, spc_basis, coeff_basis, ref_set=ref_enes)

    hf0k *= phycon.KCAL2EH

    return hf0k


def calc_channel_enes(chnl_infs, ref_ene,
                      chn_model, first_ground_model):
    """ Get the energies for several points on the reaction channel.
        The energy is determined by two different methods:
            (1) Read from the file system if chn_model == first_ground_model
            (2) Shift ene for the channel if chn_model != first_ground_model
    """

    if chn_model == first_ground_model:
        chn_enes = sum_channel_enes(chnl_infs, ref_ene, ene_lvl='ene_chnlvl')
    else:
        chn_enes1 = sum_channel_enes(chnl_infs, ref_ene, ene_lvl='ene_reflvl')
        chn_enes2 = sum_channel_enes(chnl_infs, ref_ene, ene_lvl='ene_reflvl')
        chn_enes = shift_enes(chn_enes1, chn_enes2)

    return chn_enes


def sum_channel_enes(channel_infs, ref_ene, ene_lvl='ene_chnlvl'):
    """ sum the energies
    """

    # Initialize sum ene dct
    sum_ene = {}

    # Calculate energies for species
    reac_ene = 0.0
    reac_ref_ene = 0.0
    for rct in channel_infs['reacs']:
        reac_ene += rct[ene_lvl]
        reac_ref_ene += rct['ene_tsref']
        ioprinter.info_message('reac ene', rct[ene_lvl], rct['ene_tsref'])
    sum_ene.update({'reacs': reac_ene})

    prod_ene = 0.0
    prod_ref_ene = 0.0
    for prd in channel_infs['prods']:
        prod_ene += prd[ene_lvl]
        prod_ref_ene += prd['ene_tsref']
        ioprinter.info_message('prod ene', prd[ene_lvl], prd['ene_tsref'])
    sum_ene.update({'prods': prod_ene})

    # Calculate energies for fake entrance- and exit-channel wells
    if 'fake_vdwr' in channel_infs:
        vdwr_ene = reac_ene - (1.0 * phycon.KCAL2EH)
        sum_ene.update(
            {'fake_vdwr': vdwr_ene, 'fake_vdwr_ts': reac_ene}
        )
    if 'fake_vdwp' in channel_infs:
        vdwp_ene = prod_ene - (1.0 * phycon.KCAL2EH)
        sum_ene.update(
            {'fake_vdwp': vdwp_ene, 'fake_vdwp_ts': prod_ene}
        )

    ioprinter.debug_message(
        'REAC HoF (0 K) spc lvl kcal/mol: ', reac_ene * phycon.EH2KCAL)
    ioprinter.debug_message(
        'REAC HoF (0 K) ts lvl kcal/mol: ', reac_ref_ene * phycon.EH2KCAL)
    ioprinter.debug_message(
        'PROD HoF (0 K) spc lvl kcal/mol: ', prod_ene * phycon.EH2KCAL)
    ioprinter.debug_message(
        'PROD HoF (0 K) ts lvl kcal/mol: ', prod_ref_ene * phycon.EH2KCAL)
    # Scale all of the current energies in the dict
    for spc, ene in sum_ene.items():
        sum_ene[spc] = (ene - ref_ene) * phycon.EH2KCAL

# Set the inner TS ene and scale them

    if channel_infs['ts'][0]['writer'] in ('pst_block', 'vrctst_block'):
        if len(channel_infs['reacs']) == 2:
            ts_enes = [sum(inf['ene_chnlvl'] for inf in channel_infs['reacs'])]
        else:
            ts_enes = [sum(inf['ene_chnlvl'] for inf in channel_infs['prods'])]
        channel_infs['ts'][0].update({'ene_chnlvl': ts_enes})
    else:
        if 'rpath' in channel_infs['ts']:
            ts_enes = [dct[ene_lvl] for dct in channel_infs['ts']['rpath']]
        else:
            ts_enes = [dct[ene_lvl] for dct in channel_infs['ts']]
            # ts_enes = [channel_infs['ts'][ene_lvl]]
        ioprinter.debug_message(
            'TS HoF (0 K) ts lvl kcal/mol: ', ts_enes[0] * phycon.EH2KCAL)
        if reac_ref_ene:
            if abs(ts_enes[0] - reac_ref_ene) < abs(ts_enes[0] - prod_ref_ene):
                ts_enes = [ene - reac_ref_ene + reac_ene for ene in ts_enes]
            else:
                ts_enes = [ene - prod_ref_ene + prod_ene for ene in ts_enes]
        ioprinter.debug_message(
            'TS HoF (0 K) approx spc lvl kcal/mol: ',
            ts_enes[0] * phycon.EH2KCAL)
    ts_enes = [(ene - ref_ene) * phycon.EH2KCAL for ene in ts_enes]

    sum_ene.update({'ts': ts_enes})

    return sum_ene


def shift_enes(chn_enes1, chn_enes2):
    """ When two channels dont match, the energies need to be shifted
        to bring them into alignment.
    """

    # Find a species that has enes with both methods to be used to scale
    # I don't think we need to use any species, so I will use the first
    for spc in chn_enes1:
        if chn_enes1[spc] is not None and chn_enes2[spc] is not None:
            scale_ref_spcs = spc
            break
    scale_ref_ene1 = chn_enes1[scale_ref_spcs]
    scale_ref_ene2 = chn_enes2[scale_ref_spcs]

    # Now return a dct with the lvl1 enes or the scaled lvl2 enes
    fin_enes = {}
    for spc in chn_enes1:
        if chn_enes1[spc] is not None:
            fin_enes[spc] = chn_enes1[spc]
        else:
            fin_enes[spc] = scale_ref_ene1 + (chn_enes2[spc] - scale_ref_ene2)

    return fin_enes
