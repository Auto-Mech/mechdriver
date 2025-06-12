""" Extract branching fractions of the products from
    hotenergies and pedspecies
    Fit and write the rates
"""

import sys
from scipy.interpolate import interp1d
import numpy as np
import pandas as pd
from mechanalyzer import calculator

######################## wrapper functions #########################################################


def bf_tp_dct(modeltype, ped_df, hoten_df, bf_threshold = 1e-3, rxn='', savefile=False, fne=None):
    """ Build a branching fractions dictionary as a
        function of temeprature and pressure
        containing the BFs of each product of the PES

        :param modeltype: model used for P(E1) calculations
        :type modeltype: str
        :param ped_df: dataframe[P][T] with the Series of energy distrib [en: prob(en)]
        :type ped_df: dataframe(series(float))
        :param hoten_df: hot branching fractions for hotspecies
        :type hoten_df: df[P][T]:df[allspecies][energies]}
        :param reac: label of original reactants producing hot species - only for file save
        :type hotreac: str
        :param fne: branching fractions at T,P for each product
            for the selected hotspecies
        :type fne: dataframe of series df[P][T]:series[species],
            dataframe(series(float)) (same as bf_tp_df)
        :return bf_tp_dct: branching fractions at T,P for each product
            for the selected hotspecies
        :rtype: dct{species: {pressure: (array(T), array(BF))}}
    """

    if modeltype == 'fne':
        bf_tp_df = fne
    else:
        bf_tp_df = bf_tp_df_full(ped_df, hoten_df)

    bf_tp_dct_species = bf_tp_df_todct(
        bf_tp_df, bf_threshold = bf_threshold, model=modeltype, rxn=rxn, savefile=savefile)
    _bf_tp_dct = bf_tp_dct_species

    return _bf_tp_dct


def merge_bf_ktp(bf_ktp_dct, ktp_dct, label, hotsp_dct):
    """ derive k' = bf*k and rename final ktp dictionary appropriately

        :param bf_tp_dct: branching fractions at T,P for each product
            for the selected hotspecies
        :type: dct{species: {pressure: (array(T), array(BF))}}
        :param ktp_dct: rates of the original reaction A=>B to split in A=>Pi
            with Pi = species in bf_ktp_dct (B decomposes to Pi)
        :type ktp_dct: dictionary {P: (T, k)}
        :param label: label of the original thermal reaction
        :param hotsp_dct: dictionary of hotspecies and
            corresponding fragments (if bimol)
        :type hotsp_dct: {species_unimol: [species_unimol],
                          species_bimol: [frag1, frag2], ...}
        :return rxn_ktp_dct: ktp dct of final rate constants for channels
        :rtype: {rxn: {P: (T, k)}}
    """

    ktp_dct_model_i = merge_bf_rates(bf_ktp_dct, ktp_dct)
    ktp_dct_model_i_new = rename_ktp_dct(
        ktp_dct_model_i, label, hotsp_dct)
    rxn_ktp_dct = ktp_dct_model_i_new

    return rxn_ktp_dct


def rename_ktp_dct(ktp_dct, label, hotsp_dct):
    """ rename ktp dictionary with appropriate names for prompt dissociation.
        ktp_dct.keys(): sp
        renamed_ktp_dct.keys(): rctname=>sp
        if sp is the original product, the reaction is reversible =
        :param rxn_ktp_dct: ktp dct of final rate constants for channels
        :type rxn_ktp_dct: {sp: {P: (T, k)}}
        :param label: label of the thermal reaction
        :type label: tuple ((reacs,),(prods,),(None))
        :return rename_ktp_dct: dct with new keys
        :rtype: {rxn_name: {P: (T, k)}}
    """

    renamed_ktp_dct = {}
    spc0 = list(set(hotsp_dct.keys()).intersection(label[1]))[0]
    if spc0 in ktp_dct.keys():
        renamed_ktp_dct[label] = ktp_dct[spc0]
        ktp_dct.pop(spc0)
    else:
        print('*Warning: hot spc {} not found in dct. dissociates completely?'.format(spc0))

    label1 = list(label[1])
    label1.remove(spc0)
    # other species
    for spc in ktp_dct.keys():
        frag_prods = label1 + list(hotsp_dct[spc])
        newkey = (label[0], tuple(frag_prods), label[2])
        renamed_ktp_dct[newkey] = ktp_dct[spc]

    return renamed_ktp_dct


######################## calculator functions ######################################################

def bf_tp_df_full(ped_df, hotbf_df):
    """ Build a branching fractions dataframe as a
        function of temperature and pressure containing the BFs of
        each product of the PES

        :param ped_df: dataframe[P][T] with the
            Series of energy distrib [en: prob(en)]
        :type ped_df: dataframe(series(float))
        :param hotbf_df: hot branching fractions for hotspecies
        :type hotbf_df: df[P][T]:df[allspecies][energies]}
        :return bf_tp_df: branching fractions at T,P
            for each product for the selected hotspecies
        :rtype: dataframe of series df[P][T]:series[species],
                dataframe(series(float))
    """

    # sort indexes
    ped_df = ped_df.sort_index()
    hotbf_df = hotbf_df.sort_index()
    ped_df, hotbf_df = calculator.rates.checks_temp_pressure_and_extend(
        ped_df, hotbf_df)
    # compute branching fractions
    # derive keys for dicts and complete set of T,P (should be identical)
    # temp_ped contains the desired T range; temp_hot may have a larger range
    temps, pressures = [ped_df.index, ped_df.columns]
    allspecies = hotbf_df[pressures[0]][temps[0]].columns  # extract species
    # for each T, P: compute BF
    bf_tp_df = pd.DataFrame(index=temps, columns=pressures, dtype=object)
    for temp in temps:
        for pressure in pressures:
            # extract ped and hoten by increasing index

            try:
                ped = ped_df[pressure][temp].sort_index()
            except AttributeError:
                print('empty ped at {:.0f} K and {:.1e} atm, skipping'.format(temp, pressure),)
                continue
            hoten = hotbf_df[pressure][temp].sort_index().index

            # reduce the energy range of hoten and ped
            hoten = hoten[(ped.index[0] <= hoten)*(hoten <= ped.index[-1])]
            # ped = ped[(min_en <= ped.index)*(ped.index <= max_en)]
            if len(hoten) > 3:
                ene_vect = ped.index
                ped_vect = ped.values
                # Series to allocate
                bf_series = pd.Series(0, index=allspecies, dtype=float)
                for spc in allspecies:
                    hoten_spc = hotbf_df[pressure][temp][spc][hoten]
                    f_hoten = interp1d(
                        hoten_spc.index, hoten_spc.values, bounds_error=False,
                        kind='cubic', fill_value=(hoten_spc.values[0], hoten_spc.values[-1]))
                    hoten_vect = f_hoten(ene_vect)
                    # recompute in an appropriate range
                    bf_series[spc] = np.trapz(ped_vect*hoten_vect, x=ene_vect)
                # renormalize for all species and put in dataframe
                if any(bf_series.values < 0):
                    print('Warning: found negative BFs at {:1.0f} K and {:1.1e} atm'
                          .format(temp, pressure))
                    bf_series = abs(bf_series)
                    #bf_series = abs(bf_series * np.array(bf_series > 0, dtype=int))
                bf_tp_df.at[temp, pressure] = bf_series/np.sum(bf_series.values)
                
        # if any nan: delete column
        if any(bf_tp_df.loc[temp].isnull()):
            bf_tp_df = bf_tp_df.drop(index=[temp])

    return bf_tp_df


def bf_tp_df_todct(bf_tp_df, bf_threshold = 1e-3, savefile=False, rxn='', model=''):
    """ Converts the dataframe of hot branching fractions to dictionary and
        excludes invalid BFs

        :param bf_tp_df: branching fractions at T,P for each product
            for the selected hotspecies
        :type bf_tp_df: dataframe of series df[P][T]:series[species],
            dataframe(series(float))
        :param bf_threshold: threshold to filter out all species with lower BFs
        :type bf_threshold: float
        :param rxn: name of the reaction producing hot species - only for file saving
        :type rxn: str
        :param model: type of model - only required for filename
        :type model: str
        :param savefile: save file with branching fractions?
        :type savefile: bool
        :return bf_tp_dct: branching fractions at T,P for each product
        :rtype: dct{species: {pressure: (array(T), array(BF))}} -
            same type as ktp dct
    """
    # get species
    temps, pressures = bf_tp_df.index, bf_tp_df.columns
    allspecies = bf_tp_df.iloc[0, 0].index
    bf_tp_dct_out = {}
    #set upper tol of bfrac to 1 if only 1 species available
    bftol_upper = 1 - 1e-12*(len(allspecies) > 1)
    
    # fill the dictionary:
    for spc in allspecies:
        num_data_highenough = 0
        bf_tp_dct_i = {}
        bf_df_sp_i = pd.DataFrame(
            np.zeros((len(temps), len(pressures))),
            index=temps, columns=pressures)

        for pressure in pressures:
            bf_temp = []
            temp_new = []
            for temp in temps:
                bfrac = bf_tp_df[pressure][temp][spc]
                if bftol_upper >= bfrac >= 1e-30:
                    # avoid too small values and 1 to avoid discontinuities
                    temp_new.append(temp)
                    bf_temp.append(bfrac)
                    bf_df_sp_i.at[temp, pressure] = bfrac
                # check if values is high enough
                if bfrac >= bf_threshold:
                    num_data_highenough += 1
                # condition on BF and temperature
            if bf_temp:
                bf_tp_dct_i[pressure] = (np.array(temp_new), np.array(bf_temp))

        if num_data_highenough > 0:
            bf_tp_dct_out[spc] = bf_tp_dct_i

        # write file with the BFs
        if num_data_highenough > 0 and savefile:
            bf_df_sp_i = bf_df_sp_i.reset_index()
            bf_df_sp_i = bf_df_sp_i.rename(columns = {'index': 0})
            header_label = np.array(sorted(bf_df_sp_i.columns[1:]), dtype=str)
            header_label = np.insert(header_label, 0, 'T[K]')
            labels = '\t\t'.join(header_label)
            np.savetxt(f'bf_{model}_{rxn}_{spc}.txt', bf_df_sp_i[sorted(bf_df_sp_i.columns)].values,
                       delimiter='\t', header=labels, fmt='%1.2e', comments='')

    return bf_tp_dct_out

def bf_df_fromktpdct(ktp_dct, reac_str, temps, pressures):
    """ finds reactions in ktp_dct with reactants of reac_str
        and computes the corresponding product branching fractions
        stores them in a bf dataframe

    Args:
        ktp_dct (dct): rates
        reac_str (str): species whose product branching fractions have to be found
        temps (str or numpy array): considered temperatures
        pressures (str or numpy array): considered pressures
    """

    reac_tuple = tuple(reac_str.split('+'))
    ktp_dct_ofreac = {}
    
    allspecies = []
    for key, val in ktp_dct.items():

        if key[0] == reac_tuple:
            prod = '+'.join(key[1])
            ktp_dct_ofreac[prod] = val # dct with product names as indexes
            allspecies.append(prod)

    # if allspecies is empty, print a warning
    if len(allspecies) == 0:
        print('*Warning: species {} not among reactants - branching fractions\
            asked for but not derived'.format(reac_str))
    elif len(allspecies) == 1:
        print('*Warning: species {} has only one product - {}'.format(reac_str, allspecies[0]))
        print('perhaps check if reversible rxns were ignored')
        print('if so: add thermo file')
    # turn the list of ktp_dct into series and 
    ktp_sum = dict.fromkeys(pressures)
    for pressure in pressures:
        ktot = pd.Series([0.]*len(temps), index = temps, dtype = float)
        for prod, ktp in ktp_dct_ofreac.items():
            tvect = ktp[pressure][0]
            kvect = ktp[pressure][1]
            ktseries = pd.Series(kvect, index = tvect)
            for T in tvect:
                if T not in temps:
                    ktseries = ktseries.drop(T, axis=0)
                    
            ktp_dct_ofreac[prod][pressure] = ktseries
            ktot.loc[temps] = ktot.loc[temps] + ktseries.loc[temps]
            
        ktp_sum[pressure] = ktot

    # compute branching fractions
    # for each T, P: compute BF
    bf_tp_df = pd.DataFrame(index=temps, columns=pressures, dtype=object)

    for pressure in pressures:
        for temp in temps:
            ktot_tp = ktp_sum[pressure][temp]
            bf_series = pd.Series(0, index=allspecies, dtype=float)
            for prod in allspecies:
                bf_series[prod] = ktp_dct_ofreac[prod][pressure][temp]/ktot_tp

            # warning for neg vals - absolute value
            if any(bf_series.values < 0):
                # remove warning- too verbose
                # print('Warning: found negative BFs at {:1.0f} K and {:1.1e} atm'
                #         .format(temp, pressure))
                bf_series = abs(bf_series)
            bf_tp_df.at[temp, pressure] = bf_series/np.sum(bf_series.values)
            #bf_tp_df[pressure][temp] = bf_series/np.sum(bf_series.values)

    return bf_tp_df

def merge_bf_rates(bf_tp_dct, ktp_dct):
    """ Read the branching fractions of the products and
        derive final rate constants

        :param bf_tp_dct: branching fractions at T,P for each product
            for the selected (hot)species
        :type: dct{species: {pressure: (array(T), array(BF))}}}
        :param ktp_dct: set of rates for the "total" rate constant
        :type ktp_dct: dictionary {P: (T, k)}
        :return new_ktp_dct: rates for each species multiplied by bf: bf_i*k
        :rtype: dct {species: {P: (T, k)}}
    """
    # drop "high" from ktp_dct
    if 'high' in ktp_dct.keys():
        ktp_dct.pop('high')

    # get all species and preallocate new dct
    all_species = bf_tp_dct.keys()
    new_ktp_dct = dict.fromkeys(all_species)

    # loop over all species
    for spc, bf_tp_dct_sp in bf_tp_dct.items():

        pressure_all = list(bf_tp_dct_sp.keys())
        pressure_ktp = list(ktp_dct.keys())
        # extend ktp_dct to the full range of pressures
        if not all(pressure in pressure_ktp for pressure in pressure_all):
            print('Warning: P range of ktp dictionary extended to match BF:\n')
            print(' Values at other pressures approximated from available \n')
            ktp_dct = calculator.rates.extend_dct_with_pressure(
                pressure_all, pressure_ktp, ktp_dct)

        if not all(pressure in pressure_all for pressure in pressure_ktp):
            print('Warning: P range of BF dictionary extended to match ktp:\n')
            print(' Values at other pressures approximated from available \n')
            bf_tp_dct_sp = calculator.rates.extend_dct_with_pressure(
                pressure_ktp, pressure_all, bf_tp_dct_sp)

        pressure_all = list(set(pressure_all + pressure_ktp))
        bf_ktp_dct = dict.fromkeys(pressure_all)
        for pressure in pressure_all:
            temp_bf = bf_tp_dct_sp[pressure][0]
            temp_ktp = np.array(ktp_dct[pressure][0])
            # choose the smallest vector
            temp_common = temp_bf[
                [temp_bf_i in temp_ktp for temp_bf_i in temp_bf]]
            # reconstruct appropriate datasets
            bf_series = pd.Series(bf_tp_dct_sp[pressure][1],
                                  index=temp_bf)
            ktp_series = pd.Series(np.array(ktp_dct[pressure][1]),
                                   index=temp_ktp)
            # reconstruct k and allocate
            kvals = (bf_series[temp_common].values *
                     ktp_series[temp_common].values)
            bf_ktp_dct[pressure] = (temp_common, kvals)
        # allocate the new dictionary
        new_ktp_dct[spc] = bf_ktp_dct

    return new_ktp_dct
