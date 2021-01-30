"""
  Handle labels
"""

from mechlib.amech_io import parser
from mechlib.amech_io import printer as ioprinter
from mechroutines.pf.models.typ import need_fake_wells


def make_pes_label_dct(rxn_lst, pes_idx, spc_dct, spc_model_dct):
    """ Builds a dictionary that matches the mechanism name to the labels used
        in the MESS input and output files for the whole PES
    """

    pes_label_dct = {}
    for rxn in rxn_lst:
        ioprinter.debug_message('rxn\n', rxn)
        chn_idx = rxn['chn_idx']
        pf_models = parser.model.pf_model_info(
            spc_model_dct[rxn['model'][1]]['pf'])
        ioprinter.debug_message(pf_models)
        rwell_model = pf_models['rwells']
        pwell_model = pf_models['pwells']
        tsname = 'ts_{:g}_{:g}'.format(pes_idx, chn_idx)
        pes_label_dct.update(
            _make_channel_label_dct(
                tsname, chn_idx, pes_label_dct, rxn, spc_dct,
                rwell_model, pwell_model))
        ioprinter.debug_message('pes_label dct', pes_label_dct)

    return pes_label_dct


def _make_channel_label_dct(tsname, chn_idx, label_dct, rxn, spc_dct,
                            rwell_model, pwell_model):
    """ Builds a dictionary that matches the mechanism name to the labels used
        in the MESS input and output files
    """

    # Initialize idxs for bimol, well, and fake species
    pidx, widx, fidx = 1, 1, 1
    for val in label_dct.values():
        if 'P' in val:
            pidx += 1
        elif 'W' in val:
            widx += 1
        elif 'F' in val:
            fidx += 1

    # Determine the idxs for the channel reactants
    reac_label = ''
    bimol = bool(len(rxn['reacs']) > 1)
    well_dct_key1 = '+'.join(rxn['reacs'])
    well_dct_key2 = '+'.join(rxn['reacs'][::-1])
    if well_dct_key1 not in label_dct:
        if well_dct_key2 in label_dct:
            well_dct_key1 = well_dct_key2
        else:
            if bimol:
                reac_label = 'P' + str(pidx)
                pidx += 1
                label_dct[well_dct_key1] = reac_label
            else:
                reac_label = 'W' + str(widx)
                widx += 1
                label_dct[well_dct_key1] = reac_label
    if not reac_label:
        reac_label = label_dct[well_dct_key1]

    # Determine the idxs for the channel products
    prod_label = ''
    bimol = bool(len(rxn['prods']) > 1)
    well_dct_key1 = '+'.join(rxn['prods'])
    well_dct_key2 = '+'.join(rxn['prods'][::-1])
    if well_dct_key1 not in label_dct:
        if well_dct_key2 in label_dct:
            well_dct_key1 = well_dct_key2
        else:
            if bimol:
                prod_label = 'P' + str(pidx)
                label_dct[well_dct_key1] = prod_label
            else:
                prod_label = 'W' + str(widx)
                label_dct[well_dct_key1] = prod_label
    if not prod_label:
        prod_label = label_dct[well_dct_key1]

    # Determine idxs for any fake wells if they are needed
    fake_wellr_label = ''
    if need_fake_wells(spc_dct[tsname]['class'], rwell_model):
        well_dct_key1 = 'F' + '+'.join(rxn['reacs'])
        well_dct_key2 = 'F' + '+'.join(rxn['reacs'][::-1])
        if well_dct_key1 not in label_dct:
            if well_dct_key2 in label_dct:
                well_dct_key1 = well_dct_key2
            else:
                fake_wellr_label = 'F' + str(fidx)
                fidx += 1
                label_dct[well_dct_key1] = fake_wellr_label

                # pst_r_label = 'FRB' + str(int(tsname.replace('ts_', ''))+1)
                pst_r_label = 'FRB' + str(chn_idx)
                label_dct[well_dct_key1.replace('F', 'FRB')] = pst_r_label
            if not fake_wellr_label:
                fake_wellr_label = label_dct[well_dct_key1]
                pst_r_label = label_dct[well_dct_key1.replace('F', 'FRB')]
        else:
            fake_wellr_label = label_dct[well_dct_key1]

    fake_wellp_label = ''
    if need_fake_wells(spc_dct[tsname]['class'], pwell_model):
        well_dct_key1 = 'F' + '+'.join(rxn['prods'])
        well_dct_key2 = 'F' + '+'.join(rxn['prods'][::-1])
        if well_dct_key1 not in label_dct:
            if well_dct_key2 in label_dct:
                well_dct_key1 = well_dct_key2
            else:
                fake_wellp_label = 'F' + str(fidx)
                fidx += 1
                label_dct[well_dct_key1] = fake_wellp_label

                # pst_p_label = 'FPB' + str(int(tsname.replace('ts_', ''))+1)
                pst_p_label = 'FPB' + str(chn_idx)
                label_dct[well_dct_key1.replace('F', 'FPB')] = pst_p_label
            if not fake_wellp_label:
                ioprinter.debug_message('label test', label_dct, well_dct_key1)
                fake_wellp_label = label_dct[well_dct_key1]
                if (rxn['prods'] == rxn['reacs'] or
                   rxn['prods'] == rxn['reacs'][::-1]):
                    pst_p_label = label_dct[well_dct_key1.replace('F', 'FRB')]
                else:
                    pst_p_label = label_dct[well_dct_key1.replace('F', 'FPB')]
        else:
            fake_wellp_label = label_dct[well_dct_key1]

    label_dct[tsname] = 'B' + str(chn_idx)

    return label_dct
