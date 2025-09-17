"""
    Read the product energy distribution and store the distributions obtained
"""

from doctest import OutputChecker
import sys
import numpy as np
import pandas as pd
import copy
import autoparse.find as apf
import autoparse.pattern as app
from ioformat import remove_comment_lines
from mess_io.reader._label import name_label_dct
from automol.util.dict_ import invert


def ped_names(input_str):
    """ Reads the ped_species and ped_output from the MESS input file string
        that were used in the master-equation calculation.
        :param input_str: string of lines of MESS input file
        :type input_str: str
        :return ped_species: list of names of connected REACS/PRODS
        :rtype: list(list(str))
        :return ped_output: output file with the product energy distribution
        :rtype: str
    """

    # Get the MESS input lines
    input_str = remove_comment_lines(input_str, delim_pattern=app.escape('!'))
    mess_lines = input_str.splitlines()
    check = 0
    ped_species = ()
    ped_output = None
    for line in mess_lines:
        if 'PEDSpecies' in line:
            check += 1
            ped_species_all = line.strip().split()[1:]
            for coupled_species in ped_species_all:
                newset = tuple(coupled_species.split('='))
                if newset not in ped_species:
                    ped_species += (newset,)
        if 'PEDOutput' in line:
            check += 1
            ped_output = line.strip().split()[1]
        if check == 2:
            # found everything you need: exit from loop
            break

    if not ped_species or not ped_output:
        print('*Warning: PEDSpecies and PEDOutput options incomplete, returning None \n')

    return ped_species, ped_output


def get_ped(pedoutput_str, energy_dct, sp_labels='auto'):
    """ Read `PEDOutput` file and extract product energy distribution at T,P.
        Energy in output set with respect to the ground energy of the products

        :param pedoutput_str: string of lines of ped_output file
        :type pedoutput_str: str
        :param ped_spc: species of interest in pedoutput
        :type ped_spc: list(list(str))
        :param energy_dct: energies of ped PES
        :type energy_dct: {label: energy} (str)
        :param sp_labels: type of pedspecies labels: 'inp' is how you find them
                in mess input, 'out' is how they are labeled in the output,
                'auto' sets inp if it finds the labels
        :type sp_labels: str
        :return ped_df_dct: dct(dataframe(columns:P, rows:T))
                            with the Series of energy distrib
                            for hotwells: Series of energies containing series of energy distrib for each
        :rtype ped_df_dct: {((reacs,),(prods,),(None,)): dataframe(series(float))}
        for hotwells is ((reacs,),(prods,),(None,)): dataframe(series(series((float)))
    """
    def prob_en_single(probability, energy, del_neg = False):
        # build the series and put in dataframe after
        # removing negative probs and renormalizing
        # integrate with the trapezoidal rule
        # if PEDtype is bimol, it checks for neg vals.
        prob_en = pd.Series(probability, index=energy, dtype=float)
        # remove duplicate indices 
        prob_en = prob_en[~prob_en.index.duplicated(keep = 'first')]
        prob_en = prob_en.sort_index()

        # if there are negative energies (might happen for multiple ped prods)
        prob_en = prob_en[prob_en.index > 0]

        # prob_en = prob_en[prob_en > 0]
        if del_neg:
            if np.trapz(prob_en.values, x=prob_en.index) < 0:
                # I don't know how to treat negative probabilities
                # alternative: set the absolute value? or do nothing
                prob_en *= 0
            
            elif len(prob_en[prob_en <= 0]) > 0:
                # find indexes of non positive  values
                idx_neg = np.where(prob_en <= 0)[0]
                idx_max = np.argmax(prob_en) #turn this to max val?
                # indexes of negative values closest to the center of ped (always =1)
                try:
                    low = idx_neg[idx_neg < idx_max][-1]+1
                except IndexError:
                    low = 0
                try:
                    up = idx_neg[idx_neg > idx_max][0]
                except IndexError:
                    up = -1

                prob_en = prob_en.iloc[low:up]
        #else:
        # algorithm before 05-25-2022
        # if there are negative values of the probability: remove them
        #        prob_en = prob_en[prob_en[prob_en < 0].index[-1]+1:]
        #    print(prob_en, '\n')
        # integrate with trapz
        prob_en /= abs(np.trapz(prob_en.values, x=prob_en.index))

        # if issues : keep only max value like dirac delta
        if any(np.isnan(prob_en.values)) or any(np.isinf(prob_en.values)):
            prob_en = np.nan

        try:
            if max(prob_en) > 1:
                prob_en /= max(prob_en)
        except ValueError:
            prob_max = np.argmax(probability)
            prob_en = pd.Series(probability[prob_max], index=energy[prob_max])

        except TypeError:
            pass

        return prob_en

    def indexes(label_messout, ped_lines):
        species_i = apf.where_in(label_messout, ped_lines)+1
        empty_i = apf.where_is('', ped_lines)
        final_i = np.array([empty_i[i < empty_i][0] for i in species_i])

        return species_i, final_i

    def def_prods_outinp(sp_labels, prods_list, lbl_dct):
        if sp_labels == 'inp':
            prods_outinp = pd.Series([lbl_dct[prod]
                                      for prod in prods_list], index=prods_list)
        elif sp_labels == 'out':
            prods_outinp = pd.Series(prods_list, index=prods_list)

        return prods_outinp

    def ped_and_hotwells_1(ped_lines):
        " list of ped species and hotwells for output type 1 ..."
        lines_bimolbimol, _ = indexes('Bimolecular-to-bimolecular', ped_lines)
        lines_wellbimol, _ = indexes('Well-to-bimolecular', ped_lines)
        lines_all = np.concatenate((lines_bimolbimol, lines_wellbimol))
        ped_spc = []
        for line in lines_all:
            ped_spc.extend(ped_lines[line].split()[2:])
        ped_spc = list(set(ped_spc))

        hotwells = []
        lines_wells, _ = indexes('Initial well', ped_lines)
        if len(lines_wells) > 0:
            hotwells = list(set([ped_lines[well_idx-1].split()[2]
                                 for well_idx in lines_wells]))

        return ped_spc, hotwells

    def ped_and_hotwells_2(ped_lines):
        " list of ped species and hotwells for output type 2 ..."
        lines_bimolbimol, _ = indexes('bimolecular PEDs', ped_lines)
        lines_wellbimol, _ = indexes('well PEDs', ped_lines)
        ped_spc = []
        for line in lines_bimolbimol:
            ped_spc.extend(ped_lines[line].split()[2:])
        for line in lines_wellbimol:
            ped_spc.append(
                '->'.join([ped_lines[line].split()[1], ped_lines[line+1].split()[2]]))
            # reconstruct label of the type A->B+C
        ped_spc = list(set(ped_spc))
        hotwells = []
        lines_wells, _ = indexes('hot energy', ped_lines)
        if len(lines_wells) > 0:
            hotwells = list(set([ped_lines[well_idx-1].split()[1]
                                 for well_idx in lines_wells]))
            
        return ped_spc, hotwells

    ped_lines = pedoutput_str.splitlines()

    # apf.where data of interest are
    pressure_i = apf.where_in('pressure', ped_lines)
    temperature_i = apf.where_in('temperature', ped_lines)

    # get T, P list
    pressure_lst = np.array([ped_lines[P].strip().split('=')[1]
                             for P in pressure_i], dtype=float)
    temperature_lst = np.array([ped_lines[T].strip().split('=')[1]
                                for T in temperature_i], dtype=float)
    # get label dictionary
    lbl_dct = name_label_dct(pedoutput_str)
    if lbl_dct:
        inv_lbl_dct = invert(lbl_dct)
    if sp_labels == 'auto':
        sp_labels = 'out'*(not lbl_dct) + 'inp'*(not not lbl_dct)

    outtype = 1*(len(indexes('Bimolecular-to-bimolecular', ped_lines)[0])
                 > 0) + 2*(len(indexes('bimolecular PEDs', ped_lines)[0]) > 0)

    if outtype == 1:
        ped_spc, hotwells = ped_and_hotwells_1(ped_lines)
    elif outtype == 2:
        ped_spc, hotwells = ped_and_hotwells_2(ped_lines)

    # relabel
    if sp_labels == 'inp' and len(hotwells) > 0:
        hotwells = [lbl_dct[well] for well in hotwells]

    ped_df_dct = {}

    # ped_df_dct = dict.fromkeys(energy_dct.keys())
    # find the energy for the scaling: everything refers to the products
    for label_messout in ped_spc:
        # allocate empty dataframe
        ped_df = pd.DataFrame(index=list(set(temperature_lst)),
                              columns=list(set(pressure_lst)), dtype=object)

        reacs, prods = label_messout.split('->')
        if sp_labels == 'inp':
            if 'FakeW' in lbl_dct[reacs] + lbl_dct[prods]:
                continue
            label = (tuple(lbl_dct[reacs].split('+')),
                     tuple(lbl_dct[prods].split('+')), (None,))
        elif sp_labels == 'out':
            if 'FakeW' in reacs + prods:
                continue
            label = (tuple(reacs.split('+')), 
                     tuple(prods.split('+')), (None,))
        # relabel if necessary
#        if sp_labels == 'inp':
#            label = (tuple(lbl_dct[reacs].split('+')) if 'FakeW' not in lbl_dct[reacs] else (lbl_dct[reacs],),
#                     tuple(lbl_dct[prods].split('+')) if 'FakeW' not in lbl_dct[prods] else (lbl_dct[prods],), (None,))
#        elif sp_labels == 'out':
#            label = (tuple(reacs.split('+')) if 'FakeW' not in reacs else (reacs,), 
#                     tuple(prods.split('+')) if 'FakeW' not in prods else (prods,), (None,))
        else:
            print('*Error: sp_labels must be "inp" (as in mess input) \
                or "out" (as in mess output)')
            sys.exit()

        # 0th of the energy: products energy
        prods_outinp = def_prods_outinp(sp_labels, [prods], lbl_dct)
        ene0 = -energy_dct[prods_outinp[prods]]

        species_i, final_i = indexes(label_messout, ped_lines)
        # redefine label_messout for outtype 2 for wells
        if outtype == 2 and len(label[0]) == 1:
            species_i, final_i = indexes('well: {}'.format(reacs), ped_lines) # TO BE TESTED
            species_i = np.array([i for i in species_i if 'hot' not in ped_lines[i-1]], dtype = int)+1   # remove "hot" labels
            final_i = [final_i[final_i > i][0] for i in species_i] # first empty line after each species_i
            label_messout = copy.deepcopy(prods)

        # column label
        column_i = apf.where_is(
            label_messout, ped_lines[species_i[0]-1].strip().split()[1:])[0]
        # reset species_i and final_i based on correspondence with label

        # check that length of pressure list is the same as new species_i
        # extract the data
        for i in np.arange(0, len(species_i)):
            # if labels don't match: go to next loop

            if label_messout not in ped_lines[species_i[i]-1]:
                continue

            pressure, temp = pressure_lst[i], temperature_lst[i]
            i_in, i_fin = species_i[i], final_i[i]

            en_prob_all = np.array(
                [line.strip().split() for line in ped_lines[i_in:i_fin]],
                dtype=float).T
            energy = en_prob_all[:][0] + ene0
            probability = en_prob_all[:][column_i]
            ped_df.at[temp,pressure] = prob_en_single(probability, energy)

        ped_df_dct[label] = ped_df

    for hotwell in hotwells:

        # relabel if necessary
        if sp_labels == 'inp':
            label_messout = 'well: {}'.format(inv_lbl_dct[hotwell])
        elif sp_labels == 'out':
            label_messout = 'well: {}'.format(hotwell)
        else:
            print('*Error: sp_labels must be "inp" (as in mess input) \
                or "out" (as in mess output)')
            sys.exit()
            
        # add "hot" if you have outtype 2
        if outtype == 2:
            label_messout = [label_messout, 'hot']
        # find products label
        species_i, final_i = indexes(label_messout, ped_lines)
        prods_list = ped_lines[species_i[0]].split('mol'*(outtype == 1) + 'kcal'*(outtype == 2))[1].split()

        prods_outinp = def_prods_outinp(sp_labels, prods_list, lbl_dct)

        ene0_all = pd.Series(-np.array([energy_dct[prods_outinp[prods]]
                                        for prods in prods_list]), index=prods_list)
        # allocate dataframes and labels
        for prods in prods_outinp.values:
            label = ((hotwell,), tuple(prods.split('+')), (None,))
            ped_df_dct[label] = pd.DataFrame(index=list(set(temperature_lst)),
                                             columns=list(set(pressure_lst)), dtype=object)

        # extract the data
        for i in np.arange(0, len(species_i)):
            i_in, i_fin = species_i[i]+1, final_i[i]
            pressure, temp = pressure_lst[pressure_i <
                                          i_in][-1], temperature_lst[temperature_i < i_in][-1]
        
            init_energy = float(ped_lines[i_in -
                                          2].split('=')[-1].strip().split()[0]) - energy_dct[hotwell]
            en_prob_all = np.array(
                [line.strip().split() for line in ped_lines[i_in:i_fin]],
                dtype=float).T

            for pi, prods in enumerate(prods_list):

                label = ((hotwell,), tuple(
                    prods_outinp[prods].split('+')), (None,))
                try:
                    ped_df_dct[label].at[temp, pressure][init_energy]
                except TypeError:
                    ped_df_dct[label].at[temp, pressure] = pd.Series(dtype=object)
                except KeyError:
                    pass
                finally:
                    energy = en_prob_all[:][0] + ene0_all[prods]
                    probability = en_prob_all[:][pi+1]
                    # print(probability)
                    ped_df_dct[label].at[temp, pressure][init_energy] = prob_en_single(
                        probability, energy)
                    # print(label, prods , pressure, temp , ped_df_dct[label][pressure][temp][init_energy])
                    ped_df_dct[label].at[temp, pressure] = ped_df_dct[label][pressure][temp].dropna(
                    )

    # remove "self" reactions, nonsense
    keydel = []
    for key in ped_df_dct.keys():
        if key[0] == key[1]:
            keydel.append(key)
    [ped_df_dct.pop(key) for key in keydel]
            
    return ped_df_dct
