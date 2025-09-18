"""
    Read the HotEnergies distribution and store it
    Reads mess inp and log files
"""
import sys
import numpy as np
import pandas as pd
import autoparse.find as apf
from mess_io.reader._pes import pes
from mess_io.reader._label import name_label_dct
from automol.util.dict_ import invert


def get_hot_species(input_str):
    """ Reads the HotSpecies from the MESS input file string
        that were used in the master-equation calculation.
        returns dictionary with associated energies
        :param input_str: string of lines of MESS input file
        :type input_str: str
        :return hotspecies: dictionary of hotspecies energies
        :rtype: dct{hotspecies: en}
    """

    # Get the MESS input lines
    energy_dct, _, _, _ = pes(input_str)
    mess_lines = input_str.splitlines()
    try:
        hotsp_i = apf.where_in('HotEnergies', mess_lines)[0]
        num_hotsp = int(mess_lines[hotsp_i].strip().split()[1])
        #hotspecies = [None]*num_hotsp
        hotspecies_en = {}
        for line in mess_lines[hotsp_i+1:hotsp_i+1+num_hotsp]:
            hotname = line.strip().split()[0]
            hotspecies_en[hotname] = energy_dct[hotname]
    except IndexError:
        print('*Warning: no hotspecies found, \
              returning empty dictionary \n')
        hotspecies_en = {}

    return hotspecies_en


def extract_hot_branching(hot_log_str, hotspecies_en, species_lst,
                          sp_labels='auto', filter_out01=False):
    """ Extract hot branching fractions for a single species
        :param hot_log_str: string of mess log file
        :type hot_log_str: str
        :param hotspecies_en: dct of hotspecies and corresponding energy
        :type hotspecies_en: dct{hotspecies: en}
        :param species_lst: list of all species on the PES
        :type species_lst: list
        :param sp_labels: type of species labels: 'inp' is how you find them
                in mess input, 'out' is how they are labeled in the output;
                'auto' decides 'inp' if it finds the lbl dct
        :type sp_labels: str
        :return hoten_dct: hot branching fractions for hotspecies
        :rtype hoten_dct: dct{hotspecies: df[P][T]:df[allspecies][energies]}
    """
    # get label dictionary
    lbl_dct = name_label_dct(hot_log_str)
    if lbl_dct:
        inv_lbl_dct = invert(lbl_dct)

    species_fake_dct = {sp_i: sp_i.split(
        'FakeW-')[1] for sp_i in species_lst if 'FakeW' in sp_i}

    if sp_labels == 'auto':
        sp_labels = 'inp'*(not not lbl_dct) + 'out'*(not lbl_dct)

    lines = hot_log_str.splitlines()
    # for each species: dataframe of dataframes BF[Ti][pi]
    # each of them has BF[energy][species]
    # preallocations
    hotspecies_lst = list(hotspecies_en.keys())
    # 1. extract P, T and preallocate dictionary
    pt_i_array = apf.where_in(['Pressure', 'Temperature'], lines)
    pt_list = []
    for pt_i in pt_i_array:
        pt_list.append([
            float(var)
            for var in lines[pt_i].strip().split()[2:7:4]])
    pressures = list(set([pt[0] for pt in pt_list]))
    temps = list(set([pt[1] for pt in pt_list]))
    pressures.sort(), temps.sort()

    hoten_dct = {s: pd.DataFrame(index=temps, columns=pressures)
                 for s in hotspecies_lst}

    # variables limiting the blocks
    hot_i_array = apf.where_in(['Hot distribution branching ratios'], lines)
    if len(hot_i_array) == 0:
        # different output
        hot_i_array = apf.where_in(['hot energies branching fractions'], lines)
    end_hot_i_array = apf.where_in(
        ['prompt', 'isomerization', 'dissociation'], lines)

    # 2. find Hot distribution branching ratios:
    for i, hot_i in enumerate(hot_i_array):

        # extract block, PT, and species for which BF is assigned
        lines_block = lines[hot_i+2:end_hot_i_array[i]]
        _press, _temp = [
            float(var)
            for var in lines[pt_i_array[i]].strip().split()[2:7:4]]

        # options for different outputs:
        if 'WellE' in lines[hot_i+1]:
            species_bf_i_messout = lines[hot_i+1].strip().split()[2:-1]
            outtype = 2
        elif 'kcal ' in lines[hot_i+1]:
            species_bf_i_messout = lines[hot_i+1].strip().split()[3:]
            outtype = 1
        else:
            print('*Error in reading hoten blocks - Yuri changed output again. exiting')
            sys.exit()

        # for each hotspecies: read BFs
        # rescale energy by the hotspecies energy on the PES!!
        # ref 0 energy is the ref for the PES, even for the hotspecies

        for hotspecies in hotspecies_lst:
            hot_e_lvl, branch_ratio = [], []
            ref_en = hotspecies_en[hotspecies]
            if sp_labels == 'inp':
                hotspecies_messout = inv_lbl_dct[hotspecies]
                species_bf_i = [lbl_dct[sp_i] for sp_i in species_bf_i_messout]
            elif sp_labels == 'out':
                hotspecies_messout = hotspecies
                species_bf_i = species_bf_i_messout
            else:
                raise ValueError('*Error: sp_labels must be "inp" (as in mess input) \
                    or "out" (as in mess output)')

            #species_bf_i_nofake = [sp_i for sp_i in species_bf_i if 'FakeW' not in sp_i]
            #species_bf_i_fake_dct = {}
            #for sp_i in species_bf_i:
            #    if sp_i not in species_bf_i_nofake:
            #        species_bf_i_fake_dct[sp_i] = sp_i.split('FakeW-')[1]
            
            sp_i = apf.where_is(hotspecies_messout, species_bf_i_messout)

            for line in lines_block:
                line = line.strip()
                if line.startswith(hotspecies_messout):
                    hot_e = float(line.split()[1]) - ref_en
                    # pick only values above 0
                    if hot_e not in hot_e_lvl and hot_e > 0:
                        if outtype == 1:
                            branch_ratio_arr = np.array(
                                list(line.split()[2:]), dtype=float)
                        elif outtype == 2:
                            branch_ratio_arr = np.array(
                                list(line.split()[2:-1]), dtype=float)

                        # check that value of reactant branching is between 0 and 1
                        # if any bf > 1: skip the line
                        if filter_out01:
                            if sp_i.size > 0:
                                if branch_ratio_arr[sp_i] > 1:
                                    branch_ratio_arr = np.ones(branch_ratio_arr.shape)*1e-19
                                    branch_ratio_arr[sp_i] = 1

                            # remove negative values or values >1
                            _arr = [abs(x*int(1e-20 < x <= 1))
                                    for x in branch_ratio_arr]
                            br_filter = np.array(_arr, dtype=float)

                            # if all invalid: do not save

                            if all(br_filter == 0):
                                continue
                            br_renorm = br_filter/np.sum(br_filter)
                        else:
                            br_renorm = branch_ratio_arr
                            
                        # append values                            
                        branch_ratio.append(br_renorm) ### COMMENT THIS LATER
                        hot_e_lvl.append(hot_e)

            hot_e_lvl = np.array(hot_e_lvl)
            branch_ratio = np.array(branch_ratio)

            # 3. allocate in the dataframe
            bf_hotspecies = pd.DataFrame(
                0, index=hot_e_lvl, columns=species_lst)
            bf_hotspecies[species_bf_i] = branch_ratio

            for fakesp, realsp in species_fake_dct.items():
                # sum fake species with respective wells
                bf_hotspecies[realsp] = bf_hotspecies[realsp] + bf_hotspecies[fakesp]
                bf_hotspecies = bf_hotspecies.drop(columns=[fakesp])

            hoten_dct[hotspecies].at[_temp, _press] = bf_hotspecies
 
    return hoten_dct


def extract_fne(log_str, sp_labels='auto'):
    """ Extract fne from log file
        :param log_str: string of mess log file
        :type log_str: str
        :param sp_labels: type of species labels: 'inp' is how you find them
                in mess input, 'out' is how they are labeled in the output;
                'auto' sets to inp if it finds the lbl dct, otherwise 'out'
        :type sp_labels: str
        :return dct_bf_tp_df: branching fractions at T,P
            for each product for the selected species
        :rtype: dct{sp: dataframe of series df[P][T]:series[species]},
                dct{str: dataframe(series(float))}
                made so that you can extract bf_tp_df from here to
                use it with bf_tp_df_todct
    """
    lines = log_str.splitlines()
    # get label dictionary and count N of wells
    lbl_dct = name_label_dct(log_str)
    if sp_labels == 'auto':
        sp_labels = 'inp'*(not not lbl_dct) + 'out'*(not lbl_dct)

    if sp_labels == 'inp' and lbl_dct:
        species = list(lbl_dct.values())
    elif sp_labels == 'out' and lbl_dct:
        species = list(lbl_dct.keys())
    elif sp_labels == 'out' and not lbl_dct:
        wells = [line.split('WELL: ')[1].strip()
                 for line in lines if 'WELL' in line]
        n_wells = len(wells)
        bimol = [line.split('BIMOLECULAR: ')[1].strip()
                 for line in lines if 'BIMOLECULAR' in line]
        species = wells+bimol
    else:
        raise ValueError('*Error: sp_labels must be "inp" (as in mess input) \
            or "out" (as in mess output)')

    if lbl_dct:
        n_wells = sum(
            np.array(['W' in key for key in lbl_dct.keys()], dtype=int))
        wells = species[:n_wells]

    # 1. extract P, T and preallocate dictionary
    pt_i_array = apf.where_in(['Pressure', 'Temperature'], lines)
    pt_list = []
    for pt_i in pt_i_array:
        pt_list.append([
            float(var)
            for var in lines[pt_i].strip().split()[2:7:4]])
    pressures = list(set([pt[0] for pt in pt_list]))
    temps = list(set([pt[1] for pt in pt_list]))
    pressures.sort(), temps.sort()

    # variables limiting the blocks
    fne_i_array = apf.where_in(
        ['prompt', 'isomerization', 'dissociation'], lines) + 2

    dct_bf_tp_df = dict.fromkeys(wells)

    # 2. find prompt isomerization/dissociation branching ratios:
    for i, line_in in enumerate(fne_i_array):
        line_fin = line_in + n_wells
        # extract block and PT
        lines_block = lines[line_in:line_fin]
        _press, _temp = pt_list[i]

        # for each well: read BFs
        for j, well in enumerate(wells):
            # generate df - cannot do it above otherwise overwrites
            if dct_bf_tp_df[well] is None:
                dct_bf_tp_df[well] = pd.DataFrame(
                    index=temps, columns=pressures, dtype=object)

            # read BFs
            bf_raw = lines_block[j].strip().split()[1:]
            # if pdiss=1, well values do not appear

            if len(bf_raw) == len(species) + 1:
                bf_fne = np.array(bf_raw[:n_wells] +
                                  bf_raw[n_wells+1:], dtype=float)

            elif len(bf_raw) == len(species) - n_wells + 1:
                bf_fne = np.array([0]*len(wells) +
                                  bf_raw[1:], dtype=float)
            else:
                print('I had not forseen this :C check mess log and ask Yuri')
                sys.exit()

            # remove negative values and renormalize
            if any(bf_fne < 0):
                # print removed, too verbose
                # print('Warning: found negative fne BFs at {:1.0f} K and {:1.1e} atm'
                #       .format(_temp, _press))
                # bf_fne = abs(bf_fne * np.array(bf_fne > 0, dtype=int))
                bf_fne = abs(bf_fne)
                
            bf_fne_renorm = bf_fne/np.sum(bf_fne)

            # put in dct
            dct_bf_tp_df[well].at[_temp, _press] = pd.Series(
                bf_fne_renorm, index=species, dtype=float)

    return dct_bf_tp_df
