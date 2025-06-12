""" Return the rates from a prompt dissociation process
"""
from operator import mod
import sys
import copy
import numpy
import pandas as pd
from mess_io.reader import get_species
from mess_io.reader import rates
from mess_io.reader import ped_info
from mess_io.reader import hot_info
from mechanalyzer import calculator
from mechanalyzer.calculator import thermo
from mechanalyzer.parser._util import remove_rev_rxns
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import fsolve

def prompt_dissociation_ktp_dct(ped_inp_str, ped_out_str,
                                ped_ped_str, ped_ke_out_str,
                                hot_inp_str, hot_log_str,
                                model, bf_thresh, hot_ped_str=None, hot_ke_out_str=None):
    """ Parses the MESS files and generates the rate constants
        from a prompt dissociation process
    """

    # PED INFO
    ped_spc, ped_dct, \
        dos_df, energy_dct = ped_info(
            ped_inp_str, ped_ped_str, ped_ke_out_str)
    dof_dct = calc_dof_dct(ped_inp_str, ped_spc)
    
    # HOTEN INFO
    hot_frag_dct, hot_spc_en, hoten_dct, fne_bf = \
        hot_info(hot_inp_str, hot_log_str)

    # OBTAIN ALL OF THE RATE CONSTANTS FROM THE OUTPUT FILES
    # put dictionaries together
    rxn_ktp_dct = rates.get_rxn_ktp_dct(
        ped_out_str, filter_kts=True,
        filter_reaction_types=('fake', 'self',
                                       'loss', 'capture'),
        relabel_reactions=True
    )
    # remove bw rxn based on ped_dct
    rxn_ktp_dct = remove_rev_rxns(rxn_ktp_dct, ped_dct.keys())
    # Derive Branching Fractions, Calculate Prompt Rates
    # Merge Prompt Rates with Thermal Rates
    full_prompt_rxn_ktp_dct = {}
    ene_bw_dct = {}
    pedhot_df_dct = {}

    for label, ped_df in ped_dct.items():
        modelfor = copy.deepcopy(model)
        # start prompt chains from wells/bimol with ground energies:
        # A+B=>C*+D*, (A+B=>)E=>C*+D*
        frag_reacs = label[0]
        frag_prods = label[1]
        reacs = '+'.join(frag_reacs)
        prods = '+'.join(frag_prods)
        rxn = '='.join([reacs, prods])
        ############################ various checks and derivation of frag names ###################
        # skip PEDs going to other wells or PEDs coming from hotenergy for some reason
        # - which means that the ped_df[T][p][energy] is a Series and not a float
        try:
            if len(frag_prods) < 2 or not isinstance(ped_df.iloc[0].iloc[0].iloc[0], float):
                continue
        except AttributeError:
            # something wrong with the PED
            print('*Warning: rxn {} skipped'.format(rxn))
            continue
        # select the fragment of which you want the PED:
        try:
            frag1 = list(set(hot_spc_en).intersection(frag_prods))[0]
            frag2 = list(set(frag_prods).difference((frag1,)))[0]
        except IndexError:
            print('no superposition between PED fragments {} and hot fragments {}'
                  '- skipping'.format(prods, hot_spc_en.keys()))
            continue

        ene_bw_dct[label] = energy_dct[reacs] - energy_dct[prods]
        print(
            'Processing reaction {} with DH of {:.2f} kcal/mol'.format(label, -ene_bw_dct[label]))

        if ene_bw_dct[label] < 0 and len(frag_reacs) > 1:
            print(
                'Warning: endothermic reaction {} with DH of {:.2f} kcal/mol: \
                    setting thermal model'.format(label, -ene_bw_dct[label]))
            modelfor = 'thermal'
        elif ene_bw_dct[label] < 0 and len(frag_reacs) == 1:
            print('Endothermic reaction {} with DH of {:.2f} kcal/mol excluded from loop'.format(
                label, -ene_bw_dct[label]))
            continue

        #################################################################################################

        # DERIVE PED OF THE HOT FRAGMENT
        print('deriving fragment energy distributions for {} to {}'.format(reacs, prods))
        if modelfor != 'fne':
            ped_df_frag1 = calculator.ene_partition.ped_frag1(
                ped_df, frag1, frag2, modelfor,
                dos_df=dos_df, dof_info=dof_dct[prods])
        else:
            ped_df_frag1 = None

        # calculate prompt fraction
        full_prompt_rxn_ktp_dct = calc_bf_ktp(full_prompt_rxn_ktp_dct, modelfor, bf_thresh,
                                              ped_df_frag1, fne_bf[frag1], label,
                                              hoten_dct[frag1], rxn, rxn_ktp_dct[label], hot_frag_dct)
        # print(full_prompt_rxn_ktp_dct.keys(), '\n')
        # if hot options: it means hoten has also ped - extract info

        if hot_ke_out_str and hot_ped_str:
            # frag1 hotspecies - devi prendere quella distribuzione dalle hot
            pedhot_df_dct_spc, ene_bw_spc = build_pedhot_df_dct(hot_inp_str, hot_ped_str, hot_ke_out_str,
                                                                frag1, ped_df_frag1, ene_bw_dct[label], model, label)
            pedhot_df_dct.update(pedhot_df_dct_spc)
            ene_bw_dct.update(ene_bw_spc)

    return full_prompt_rxn_ktp_dct, pedhot_df_dct, ene_bw_dct


def prompt_chain_ktp_dct(rxn_ktp_dct, ene_start_df_dct,
                         pedhot_inp_str, pedhot_log_str,
                         model, bf_thresh, ene_bw_dct, hot_ped_str=None,
                         hot_ke_out_str=None):
    """ Starts from rxn_ktp_dct
        rxn_ktp_dct: ktp dct to update with prompt fractions
        ene_start_df_dct: starting energy distribution for given reaction chain
                     {(reacs,),(prods,),(None,): {frag1: energy_distr, frag2: energy_distr}

    """
    # HOTEN INFO
    hot_frag_dct, hot_spc_en, hoten_dct, fne_bf = \
        hot_info(pedhot_inp_str, pedhot_log_str)

    # Derive Branching Fractions, Calculate Prompt Rates
    # Merge Prompt Rates with Thermal Rates
    full_prompt_rxn_ktp_dct = {}
    # print(ene_start_df_dct.keys(), hot_spc_en.keys())
    for label, ene_start_df in ene_start_df_dct.items():
        modelfor = copy.deepcopy(model)
        # find dissociating fragment if available
        frag10, frag20 = ene_start_df.keys()
        frag = frag10*(frag10 in hot_spc_en.keys()) + \
            frag20*(frag20 in hot_spc_en.keys())

        # continue if no overlap
        if not frag:
            continue

        rxn = '+'.join(label[0])+'='+'+'.join(label[1])
        print(
            'Processing reaction {} with DH of {:.2f} kcal/mol'.format(label, -ene_bw_dct[label]))
        # model checks
        if ene_bw_dct[label] < 0:
            print(
                'Warning: endothermic reaction {} with DH of {:.2f} kcal/mol: \
                   setting thermal model'.format(label, -ene_bw_dct[label]))
            modelfor = 'thermal'

        full_prompt_rxn_ktp_dct = calc_bf_ktp(full_prompt_rxn_ktp_dct, modelfor, bf_thresh,
                                              ene_start_df[frag], fne_bf[frag], label,
                                              hoten_dct[frag], rxn, rxn_ktp_dct[label], hot_frag_dct)

        # if hot options: it means hoten has also ped - extract info
        pedhot_df_dct = {}
        if hot_ke_out_str and hot_ped_str:
            # frag1 è hotspecies - devi prendere quella distribuzione dalle hot
            pedhot_df_dct_spc, ene_bw_spc = build_pedhot_df_dct(pedhot_inp_str, hot_ped_str, hot_ke_out_str,
                                                                frag, ene_start_df[frag], ene_bw_dct[label], model, label)
            pedhot_df_dct.update(pedhot_df_dct_spc)
            ene_bw_dct.update(ene_bw_spc)

    return full_prompt_rxn_ktp_dct, pedhot_df_dct, ene_bw_dct


def build_pedhot_df_dct(hot_inp_str, hot_ped_str, hot_ke_out_str,
                        starthotfrag, starthotfrag_df, ene_bw, model, label_start):
    """ label_start: label of the starting reaction (reaction generating the radical)
    """

    pedhot_df_dct_tot = {}
    pedhot_df_dct = {}
    ene_bw_dct = {}
    hot_ped_spc, hot_ped_dct, \
        hot_dos_df, hot_energy_dct = ped_info(
            hot_inp_str, hot_ped_str, hot_ke_out_str)
    hot_dof_dct = calc_dof_dct(hot_inp_str, hot_ped_spc)
    # starting energy distribution
    starthot_df = starthotfrag_df

    # just 1 ped set in this case
    _, prods = hot_ped_spc[0]
    frag1, frag2 = prods.split('+')
    label = ((starthotfrag,), tuple(prods.split('+')), (None,))

    newprods = tuple([pr for pr in label_start[1] if pr !=
                      starthotfrag]) + tuple(prods.split('+'))

    newlabel = (label_start[0], newprods, label_start[2])
    ene_bw += (hot_energy_dct[starthotfrag] - hot_energy_dct[prods])

    print('deriving hot fragment distribution for {}'.format(newlabel))
    if ene_bw < 0:
        print(
            'endothermic reaction {} with DH of {:.2f} kcal/mol: \
                setting thermal model'.format(newlabel, -ene_bw))
        model = 'thermal'

    # calculate new PED (beta)
    ped_df_fromhot = hot_ped_dct[label]

    ped_df_rescaled = calculator.ene_partition.ped_df_rescale(
        starthot_df, ped_df_fromhot, save = True, name = '+'.join(newlabel[0]) + '=' + '+'.join(newlabel[1]))

    ########################################################################################
    # stupid test to delete later: try to just have the starthot_df but rescale the energy
    """
    ped_df_rescaled = calculator.ene_partition.ped_df_rescale_test(
        starthot_df, hot_energy_dct[starthotfrag]-hot_energy_dct[prods], save=boolsave)
    """
    #########################################################################################

    # assign name to the global reaction
    if model != 'fne':

        # DERIVE PED OF THE HOT FRAGMENT - BOTH , CHECK FRAG2 NOT ATOM
        pedhot_df_dct[frag1] = calculator.ene_partition.ped_frag1(
            ped_df_rescaled, frag1, frag2, model,
            dos_df=hot_dos_df, dof_info=hot_dof_dct[prods])

        if hot_dof_dct[prods]['n_atoms'][frag2] > 1:
            pedhot_df_dct[frag2] = calculator.ene_partition.ped_frag1(
                ped_df_rescaled, frag2, frag1, model,
                dos_df=hot_dos_df, dof_info=hot_dof_dct[prods])

    pedhot_df_dct_tot[newlabel] = copy.deepcopy(pedhot_df_dct)
    ene_bw_dct[newlabel] = ene_bw

    return pedhot_df_dct_tot, ene_bw_dct


def calc_bf_ktp(full_prompt_rxn_ktp_dct, model, bf_thresh,
                ped_df_frag1, fne_bf_frag1, label,
                hoten_dct_frag, rxn, ktp_dct, hot_frag_dct):
    """ generate full bf and ktp dct from ped and hoten
        NB frag2 must be tuple
    """
    # JOIN PED AND HOTEN -> DERIVE PRODUCTS BF
    bf_tp_dct = calculator.bf.bf_tp_dct(
        model, ped_df_frag1, hoten_dct_frag, bf_threshold = bf_thresh,
        savefile=True, rxn=rxn, fne=fne_bf_frag1)

    # CALCULATE PROMPT DISSOCIATION RATES K*BF

    prompt_rxn_ktp_dct = calculator.bf.merge_bf_ktp(
        bf_tp_dct, ktp_dct,
        label, hot_frag_dct)

    # Merge Prompt Rates with all current; Rates added if rxn is prev. found
    # print(modeltype, full_prompt_rxn_ktp_dct[modeltype], '\n', prompt_rxn_ktp_dct[modeltype], '\n','\n','stop')
    full_prompt_rxn_ktp_dct = calculator.rates.merge_rxn_ktp_dcts(
        full_prompt_rxn_ktp_dct,
        prompt_rxn_ktp_dct
    )

    return full_prompt_rxn_ktp_dct

def calc_dof_dct(ped_inp_str, ped_spc):
    spc_blocks_ped = get_species(ped_inp_str)
    dof_dct = {}
    for spc in ped_spc:
        _, prods = spc

        # Derive dofs involved
        dof_dct[prods] = calculator.spinfo_frommess.get_info(
            spc_blocks_ped[prods])
        
    return dof_dct


###################### FUNCTIONS FOR THE SORTER #######################
################### these work with dataframes ########################
def get_max_reactivity(hot_sp, hot_sp_df, therm_df, T0, Tref):
    """ input: dataframe of decomposition reactions of a species
        return the maximum deco rate k[Tref] and the minimum dh[T0] for dissociation
    """

    dhmin = 1e8  # put unphysical value for sure higher than all
    kmax = -1e-6
    allunimol = all([len(rcts) == 1 and len(prds) ==1 for rcts in hot_sp_df['rct_names_lst'].values 
                     for prds in hot_sp_df['prd_names_lst'].values])

    for rxn in hot_sp_df.index:
        rcts = hot_sp_df['rct_names_lst'][rxn]
        prds = hot_sp_df['prd_names_lst'][rxn]
        # used to filter out too fast isomerization channels, but causes failure when only unimol channels are present
        if not allunimol and ((hot_sp in prds and (len(rcts) != 2 or len(prds) != 1)) 
            or (hot_sp in rcts and (len(prds) != 2 or len(rcts) != 1))):
            # consider only unimol chnls to bimol species
            continue
        
        label = rxn[0]
        
        # analyze DH
        dh_rxn = thermo.extract_deltaX_therm(therm_df, rcts, prds, 'H')/1000
        if hot_sp in prds:  # revert sign
            dh_rxn = -dh_rxn
            label = '{}={}'.format(rxn[0].split('=')[-1], rxn[0].split('=')[0])
            
        if dh_rxn[T0] < dhmin:
            dh_min_hot = copy.deepcopy(dh_rxn)
            dhmin = copy.deepcopy(dh_rxn[T0])
        
        # analyze k
        k_series = pd.Series(
            hot_sp_df['param_vals'][rxn][0][1.][1], index=hot_sp_df['param_vals'][rxn][0][1.][0])
        if Tref not in k_series.index:
            print('*Warning: {} K not found - k not extracted'.format(Tref))
            k = numpy.NaN
            continue
        else:
            k = k_series[Tref]

        if hot_sp in prds:  # get bw rate constant
            dg_rxn = thermo.extract_deltaX_therm(
                therm_df, rcts, prds, 'G')
            # get common T values of dg and k
            k_series = calculator.rates.get_bw_rate(k_series, rcts, prds, dg_rxn)
            k = k_series[Tref]

        if k > kmax:
            k_max_hot = copy.deepcopy(k_series)
            kmax = copy.deepcopy(k)
            rxn_label = copy.deepcopy(label)

    return k_max_hot, dh_min_hot, rxn_label

def estimate_hot_hk(dh_pi, Tref, cp_hot, kmax_hot, dhmin_hot):
    """ estimate dhtot and k* of hot products for non boltzmann PES
        units cal, mol, K
        NB dh_pi is the dh transferred to the hot species
    """
    # approximate with close Tref if absent
    if dh_pi[Tref] < 0:
        T_star, k_star = kt_star(-dh_pi[Tref], cp_hot, Tref, kmax_hot)
    else:
        T_star = Tref
        k_star = kmax_hot[Tref]

    dh_tot = (dh_pi + dhmin_hot)/1000
    return T_star, k_star, dh_tot

def kt_star(DH0, Cp, T0, k_series):
    """
    Take a DH0 and Cp(T) and find the T* where DH = sum(Cp*dT) from T0 to T*
    Then compute k* at the new T*
    """
    def tozero(T_star, DH0, f_cp, T0):
        int_fct, _ = quad(f_cp, T0, T_star)
        return DH0 - int_fct
    
    try:
        f_cp = interp1d(numpy.array(Cp.index, dtype=float), numpy.array(Cp.values, dtype=float), kind='cubic')
        f_k = interp1d(numpy.array(k_series.index, dtype=float), numpy.array(k_series.values, dtype=float), kind='cubic')
    except TypeError as e:
        print('error in interpolation for Cp: {}'.format(e))
        print('Setting T* as T0 : {} K'.format(T0))
        return T0, k_series[T0]
    except ValueError as e:
        print('Interpolation failed - Cp vector must have at least four values. T vals are {}'.format(list(Cp.index)))
        print('Setting T* as T0 : {} K'.format(T0))
        print(e)
        return T0, k_series[T0]
    
    try:
        T_star = fsolve(tozero, T0+100, args=(DH0, f_cp, T0))[0]
    except ValueError:
        print(
            '*Warning: fsolve failed to compute T*, out of range - using fixed Cp at {:.0f} K'.format(T0))
        T_star = T0+DH0/Cp[T0]
    try:
        k_star = f_k(T_star)
    except ValueError:
        k_series = k_series.sort_index()
        k_star = k_series.iloc[-1]
        print(
            '*Warning: k* out of range at T of {} K - determine k at max T of {} K'.format(T_star, k_series.index[-1]))

    return T_star, k_star

