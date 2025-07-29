"""
Sorter module - sorting of the mechanism according to various options
- reaction classes
- pes/subpes
- multiplicity
- species subsets
- submechanism
"""


import sys
import copy
import pandas as pd
import numpy
import automol
from mechanalyzer.builder import submech
from mechanalyzer.builder import rxnclass
from mechanalyzer.builder import connect_rxn_df
from mechanalyzer.builder import add_wellskip
from mechanalyzer.calculator import rates as calc_rates
from mechanalyzer.calculator import thermo
from mechanalyzer.calculator.ene_partition import phi_equip_fromdct
from mechanalyzer.calculator import nonboltz
from mechanalyzer.parser import pes
from mechanalyzer.parser.spc import name_inchi_dct
from mechanalyzer.parser._util import count_atoms
from mechanalyzer.parser._util import order_rct_bystoich
from mechanalyzer.parser._util import extract_spc
from mechanalyzer.parser._util import get_mult
from mechanalyzer.parser._util import get_fml

# species groups assigned when processing submech
FUELGROUP = ['FUEL','FUEL_RAD','FUEL_ADD_H','FUEL_ADD_CH3','FUEL_ADD_O','FUEL_ADD_OH','FUEL_ADD_O2','R_CH3','R_O','R_O2','R_O4','R_O3-H']
SUBMECHGROUP = ['CORE']
SUPMECHGROUP = ['SUPFUEL','SUBFUEL']

def mech_info(rxn_param_dct, spc_dct):
    """ Build mech_info object for mech sorting

        :param spc_dct: species dictionary
        :type spc_dct: dict[?:?]
        :param rxn_dct: parameter dictionary
        :type rxn_dct: dict[?:?]
        :return mech_info: objects with mech info
        :rtype: list
    """

    def _check_names(rct_names, prd_names, thrdbdy_lst, all_spc_names):
        """ Assess if the reactant and product names provided in the
            rxn_param_dct exist in the spc_dct
        """
        # extract third body names
        # [('+M',), ('(+HE)',), ('(+M)',), (None,), (None,), ('(+M)',)]
        thrdbdy_names = ()
        for thrd in thrdbdy_lst:
            for el in thrd:
                if el is not None:
                    thrdbd = el.split('+')[-1].strip().split(')')[0].strip()
                    if thrdbd != 'M' and thrdbd not in thrdbdy_names:
                        thrdbdy_names += (thrdbd, )

        all_mech_names = ()
        for _rct_names, _prd_names in zip(rct_names, prd_names):
            all_mech_names += _rct_names
            all_mech_names += _prd_names
        all_mech_names = set(all_mech_names + thrdbdy_names)

        missing_names = all_mech_names - all_spc_names
        if missing_names:
            print('Names in provided in mechanism, '
                  'but not provided in species list (likely from .csv file):')
            for name in missing_names:
                print('  ', name)
            print('Unable to finish parsing mechanism. Exiting...')
            sys.exit()

        return all_mech_names

    def _inf(rct_names, prd_names, ich_dct):
        """ Sort reactant and product name lists by formula to facilitate
            multichannel, multiwell rate evaluations
        """
        rxn_name_lst, formula_str_lst, formula_dct_lst = [], [], []
        for _rct_names, _prd_names in zip(rct_names, prd_names):
            rxn_name = '='.join(['+'.join(_rct_names), '+'.join(_prd_names)])
            rxn_name_lst.append(rxn_name)
            rct_ichs = list(map(ich_dct.__getitem__, _rct_names))
            formula_dct, formula_str = get_fml(rct_ichs)
            formula_dct_lst.append(formula_dct)
            formula_str_lst.append(formula_str)

        return formula_dct_lst, formula_str_lst, rxn_name_lst

    if not all([isinstance(val, dict) or isinstance(val, list) for val in rxn_param_dct.values()]):
        print(
            '*ktp dct vals not found - derived for sorting purposes derived at [300:10:2010] K at 1 atm')
        rxn_ktp_dct = calc_rates.eval_rxn_param_dct(rxn_param_dct, [numpy.arange(300, 2010, 10)],
                                                    [1])
        for key, val in rxn_ktp_dct.items():
            if isinstance(val, dict):
                rxn_ktp_dct[key] = [val]
    else:
        rxn_ktp_dct = rxn_param_dct  # it means you already provided a ktp dct as input

    # Extract info from dictionary
    rcts, prds, thrdbdy = zip(*rxn_param_dct.keys())
    rct_names, prd_names, thrdbdy_lst = list(rcts), list(prds), list(thrdbdy)
    # Check if the rxn_param dct may be fully parsed
    all_mech_names =_check_names(rct_names, prd_names, thrdbdy_lst, set(spc_dct.keys()))
    # filter species dictionary to make processing lighter from now on
    spc_dct ={k: spc_dct[k] for k in all_mech_names}

    # formulas and reaction names (repplace with the mech info from ckin
    ich_dct = name_inchi_dct(spc_dct)
    formula_dct, formula_str, rxn_name = _inf(rct_names, prd_names, ich_dct)

    return [spc_dct, formula_dct, formula_str,
            rct_names, prd_names, thrdbdy_lst,
            rxn_name, list(rxn_ktp_dct.values())]


def cmts_string(name, label, cltype):
    """ assign comment strings based on reaction class

    :param name: reaction class (list, string, int, float)
    :param label: labels corresponding to reaction class (list)
    :param cltype: format type. options: class_head, class, subclass (string)

    :returns: comments string for writing in the mechanism
    :rtype: str
    """
    # assing top headers:
    head = '!!!!!!!!!!!!!!!!!!!!!!!!!\n'

    try:
        labeldct = dict(zip(name, label.values))
    except TypeError:
        labeldct = dict(zip([name], [label.values]))

    def formatname(name):
        if isinstance(name, int):
            name = str(name)
        elif isinstance(name, numpy.int64):
            name = str(name)
        elif isinstance(name, float):
            if labeldct[name] not in ['maxval', 'maxratio']:
                name = f'{name:.2f}'
            else:
                name = f'{name:.2e}'

        return name

    if isinstance(name, tuple):
        namenew = []
        for _nm in name:
            namenew.append(formatname(_nm))
        name = namenew
    else:
        name = [formatname(name)]

    if cltype == 'class_head':
        cmtlabel = '!       '+'.'.join(label)
        cmtlabel += '!       '+'.'.join(name)+'\n'
        rxnclass = head + cmtlabel + head
    else:
        cmtlabel = '.'.join(label) + '  '
        rxnclass = '! ' + cmtlabel + '.'.join(name)

    return rxnclass


def conn_chn_df(mech_df): # moved outside of the class - taken for granted: classification always includes pes, subpes, chnl
    """ Identifies connected channels and assigns them to the same subpes
        Generate pes dictionary for each reaction and save for later use

    :param self.mech_df: dataframe with mech info (contains all reactions)
    :param conn_chn_df: empty dataframe df[subpes][rxn]

    :returns: conn_chn_df dataframe[subpes][rxn]
    :rtype: dataframe[int][tuple]
    """
    conn_chn_df = pd.DataFrame(
            index = mech_df.index, columns=['subpes'])

    pes_dct_df = pd.DataFrame(
        index=mech_df.index,
        columns=['chnl', 'pes_chnl_tuple'],
        # columns=['pes_dct', 'pes_chnl_tuple', 'pes_chnl'],
        dtype=object)

    for _, peslist in mech_df.groupby('pes'):
        idx_start = 0
        # Set the names lists for the rxns and species needed below
        # print(peslist)
        peslist = peslist.sort_values(by=['rxn_names'])
        pes_rct_names_lst = peslist['rct_names_lst'].values
        pes_prd_names_lst = peslist['prd_names_lst'].values
        pes_rxn_name_lst = peslist.index
        connchnls = pes.find_conn_chnls(
            pes_rct_names_lst, pes_prd_names_lst, pes_rxn_name_lst)
        # print(connchnls)
        # Write subpes in conn_chn_df
        for key, value in connchnls.items():
            rxns = peslist.iloc[value].sort_values(
                by=['rxn_names'], ascending=False).index

            # conn_chn_df['subpes'][rxns] = key+1
            conn_chn_df.loc[rxns, 'subpes'] = key+1
            # reorder by rxn name before assigning the channel index

            for chnl_idx, rxn in enumerate(rxns):
                rxn_name = tuple(
                    (peslist['rct_names_lst'][rxn],
                        peslist['prd_names_lst'][rxn]))
                pes_dct_df.at[rxn, 'chnl'] = chnl_idx+idx_start+1

                pes_dct_df.at[rxn, 'pes_chnl_tuple'] = (
                    chnl_idx+idx_start, rxn_name)

            idx_start += len(rxns)

    conn_chn_df = pd.concat([conn_chn_df, pes_dct_df], axis=1)

    return conn_chn_df

class SortMech:
    """ class of methods to organize the mechanism according to given criteria
    """

    def __init__(self, rxn_param_dct, spc_dct):
        """ Initializes the mechanism dataframe and the species dictionary

        :param rxn_param_dct: rxn param dct for info extraction
        :param spc_dct: species dictionary

        :returns: None, updates self.
                    self.mech_df: dataframe with mech info
                    self.spc_dct: species dictionary
        """

        # Extract data from mech info
        [spc_dct, formula_dct_lst, formulas, rct_names_lst,
            prd_names_lst, thrdbdy_lst, rxn_name_lst, param_vals] = mech_info(rxn_param_dct, spc_dct)

        rxn_index = list(zip(rxn_name_lst, thrdbdy_lst))

        # Set dataframe: extract useful info
        pes_lst = count_atoms(formula_dct_lst)
        isthrdbdy = numpy.array(
            [(t[0] is not None and '(' not in t[0]) for t in thrdbdy_lst],
            dtype=int)
        molecularity = numpy.array(
            list(map(len, rct_names_lst)), dtype=int) + isthrdbdy
        n_of_prods = list(map(len, prd_names_lst))
        rct_names_lst_ordered = order_rct_bystoich(
            rct_names_lst, spc_dct=spc_dct)  # put heavier reactant first
        prd_names_lst_ordered = order_rct_bystoich(
            prd_names_lst, spc_dct=spc_dct)  # put heavier product first
        rct_1, rct_2 = extract_spc(rct_names_lst_ordered)
        data = numpy.array(
            [rct_names_lst, prd_names_lst,
             rct_names_lst_ordered, prd_names_lst_ordered,
             rct_1, rct_2, molecularity, n_of_prods, pes_lst, formulas,
             isthrdbdy, thrdbdy_lst, param_vals, rxn_name_lst],
            dtype=object).T
        self.mech_df = pd.DataFrame(
            data, index=rxn_index,
            columns=['rct_names_lst', 'prd_names_lst', 'rct_names_lst_ord',
                     'prd_names_lst_ord', 'r1', 'r2', 'molecularity',
                     'N_of_prods', 'pes_fml', 'formulas', 'isthrdbdy',
                     'thrdbdy', 'param_vals', 'rxn_names'], dtype=object)

        # reindex pes
        df_pes = pd.DataFrame(
            index=self.mech_df.index, columns=['pes'])
        pes_index = 1
        for _, pes_fml in self.mech_df.groupby('pes_fml'):
            rxns = pes_fml.index
            df_pes.loc[rxns,'pes'] = pes_index
            pes_index += 1
        # Concatenate the new portion of dataframe
        self.mech_df = pd.concat([self.mech_df, df_pes], axis=1)
        # add subpes and chnl once and for all
        self.mech_df = pd.concat(
            [self.mech_df, conn_chn_df(self.mech_df)], axis=1)  # add subpes

        self.spc_dct = spc_dct  # set for later use
        # empty list for initialization (otherwise pylint warning)
        self.species_subset_df = ()
        self.species_list = ()

    def sort(self, hierarchy, species_list):
        """ Main flow of the sorter: takes a set of reactions and classifies
            them as indicated in hierarchy; possibly filters the species
            according to the species list

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.spc_dct: species dictionary
        :param hierarchy: list of hierarchical criteria for mech organization
            ['..','..','N']: N is the N of criteria to use for headers
            all the other criteria are written as inline comments
        :param species_list: list of species subset; may be empty
            if ['onename'] and 'submech' in criteria: extracts fuel submechanism
        :returns: None. updates self
        """
        # set labels for all the possible criteria
        # maybe turn into dataframe with criteria, labels, fct names, and ascending opts
        criteria_all = ['molecularity', 'N_of_prods', 'species', 'fml',
                        'pes', 'subpes', 'chnl',
                        'r1', 'mult',
                        'rxn_class_broad', 'rxn_class_graph',
                        'submech', 'submech_prompt', 'submech_keepsubfuel',
                        'submech_deletelarge',  'rxn_max_vals', 'rxn_max_ratio',
                        ]
        labels_all = ['NR', 'N_of_prods', 'SPECIES', 'formulas',
                      'pes', 'subpes', 'channel',
                      'Heavier_rct', 'Totalmultiplicity',
                      'rxntype', 'rxnclass',
                      'submech', 'submech_prompt', 'submech_keepsubfuel',
                      'submech_deletelarge', 'maxval', 'maxratio',
                      ]
        # check that all criteria in hierarchy are allowed
        notavail = [hr for hr in hierarchy[:-1] if hr not in criteria_all]
        assert len(notavail) == 0, (
            '*Error: classification according to {} not available; criteria are: \n {}'.format(notavail,
                                                                                               '\n'.join(criteria_all))
        )
        # check on pes/subpes criteria:
        # if subpes alone, also add pes
        if (('subpes' in hierarchy and 'pes' not in hierarchy) or
                ('chnl' in hierarchy and 'pes' not in hierarchy and
                 'subpes' not in hierarchy)):
            try:
                idx_insert_pes = hierarchy.index('subpes')
            except ValueError:
                idx_insert_pes = hierarchy.index('chnl')
            hierarchy.insert(idx_insert_pes, 'pes')
            # chenge N of headers if necessary
            hierarchy[-1] += 1*(hierarchy[-1] > idx_insert_pes)
            print('pes criterion added: necessary for subpes/chnls')

        elif 'submech_prompt' in hierarchy:
            # add pes, subpes, chnl:
            hierarchy_new = ['pes', 'subpes', 'chnl']
            [hierarchy_new.append(h)
                for h in hierarchy if h not in hierarchy_new]
            hierarchy = copy.deepcopy(hierarchy_new)
            print('final hierarchy for classification: {}'.format('-'.join(hierarchy[:-1])))

        # now that you fixed it, save it in self
        self.hierarchy = hierarchy

        # if species_list is not empty: pre-process the mechanism
        self.preproc_specieslist(species_list)

        # sort and label
        self.sort_and_label(criteria_all, labels_all)

        self.criteria_all = criteria_all
        self.labels_all = labels_all

    def sort_and_label(self, criteria_all, labels_all):
        # 0. Look for keywords and sort accordingly
        # List of available sorting options with specific sorting functions
        sort_optns_dct = {
            'species': self.group_species,
            'mult': self.reac_mult,
            'rxn_class_broad': self.rxnclass_broad,
            'rxn_class_graph': self.rxnclass_graph,
            'submech': self.group_submech,
            'submech_prompt': self.group_prompt,
            'submech_keepsubfuel': self.group_submech,
            'submech_deletelarge': self.group_submech,
            'rxn_max_vals': self.rxn_max_vals,
            'rxn_max_ratio': self.rxn_max_ratio
        }

        for optn, fun_name in sort_optns_dct.items():
            # Call sorting function
            if any(optn == inp_crt for inp_crt in self.hierarchy):
                # Generate corresponding dataframe
                df_optn = pd.DataFrame(
                    index=self.mech_df.index, columns=[optn])
                df_optn = fun_name(df_optn)
                # Concatenate the new portion of dataframe
                self.mech_df = pd.concat([self.mech_df, df_optn], axis=1)

        # 0. remove fake rxns added through wellskipping generator
        rxns_fake = self.mech_df[self.mech_df['chnl']
                                == 'WELLSKIPPING FAKE'].index
        self.mech_df = self.mech_df.drop(index=rxns_fake)

        # 1. Sort
        # series: ascending/descending values
        asc_val = [True]*len(criteria_all + ['rxn_names'])
        # rxn vals, ratio and names are descending
        asc_val[-3:] = [False, False, False]
        asc_series = pd.Series(asc_val, index=criteria_all + ['rxn_names'])


        try:
            # last "standard" criterion: rxn name
            asclst = list(
                asc_series[self.hierarchy[:-1] + ['rxn_names']].values)
            self.mech_df = self.mech_df.sort_values(
                by=self.hierarchy[:-1] + ['rxn_names'],
                ascending=asclst)
        except KeyError as err:
            raise KeyError(
                'Error: Reactions not sorted according ',
                f'to all criteria: missing {err}, exiting') from err

        # 2. assign class headers
        labels = pd.Series(labels_all, index=criteria_all)
        self.class_headers(self.hierarchy, labels)

    def preproc_specieslist(self, species_list):

        sumbech_optns_dct = {'submech': {'fun_name': submech.species_subset_fuel,
                                         'filtertype': 'submech'},
                             'submech_keepsubfuel': {'fun_name': submech.species_subset_keep,
                                             'filtertype': 'submech_keepsubfuel'},
                             'submech_deletelarge': {'fun_name': submech.species_subset_del,
                                             'filtertype': 'submech_deletelarge'},
                             'submech_prompt': {'filtertype': 'submech_prompt'}}

        submech_name = [optn for optn in self.hierarchy[:-1]
                            if 'submech' in optn]
        if len(submech_name) > 0:
            submech_name = submech_name[0]

        # extract filters 'deleteabove', 'keepbelow', 'singlespecies
        stoich = None
        species_list_new = []
        for sp in species_list:
            if 'deleteabove' in sp:
                stoich = sp.split()[1]
                submech_name = 'submech_deletelarge' # set in case you don't have it in the criteria
            elif 'keepbelow' in sp:
                stoich = sp.split()[1]
                submech_name = 'submech_keepsubfuel' # set in case you don't have it in the criteria
            elif 'singlespecies' in sp:
                submech_name = 'singlespecies'
            else:
                species_list_new.append(sp)

        species_list = species_list_new

        # submech if 1 species selected and no submech_name /  = True found: select submech anyway
        if len(submech_name) == 0 and len(species_list) == 1:
            print('Species selected but no sorting criterion applied - submech is applied by default')
            print('If you want only reactions where {} appears, add "singlespecies = True "'.format(species_list[0]))
            print(' in sort_mech section "singlespecies" in the sort_lst input')
            submech_name = 'submech'

        elif len(species_list) > 1 and len(submech_name) == 0:
            print('Multiple species selected and no submech specifications: will extract reactions where {} are involved'.format(species_list))
            submech_name = 'singlespecies'

        elif len(species_list) > 1 and submech_name != 'submech_prompt':
            raise ValueError('Error: multiple species only available for submech_prompt sorting')

        elif len(species_list) == 0 and submech_name == 'submech_prompt':
            # species list includes all radicals in the mech
            print('WARNING: Prompt selected w/o species specification: \
                all radicals analyzed ...')
            for sp_i in self.spc_dct.keys():
                mult = get_mult(sp_i, self.spc_dct)
                atoms = sum(automol.chi.formula(
                    self.spc_dct[sp_i]['inchi']).values())
                if mult > 1 and atoms > 2:
                    species_list.append(sp_i)


        if submech_name in ['submech_deletelarge', 'submech_keepsubfuel']:
            # Select subset of species according to stoichiometries
            # specified in submech.py
            fuel = species_list[0] if len(species_list) == 1 else None
            print('{} extracted according to input stoich specification {} and fuel {}'.format(submech_name, stoich, fuel))
            species_list, species_subset_df = sumbech_optns_dct[submech_name]['fun_name'](
                self.spc_dct, fuel = fuel, stoich = stoich)
            self.species_subset_df = species_subset_df

        elif submech_name == 'submech' and len(species_list) == 1:
            # 'submech' option (either declared explicitly or implicitly by indicating one species to isolate)
            print('Submech extracted for fuel {}'.format(species_list[0]))
            species_list, species_subset_df, _ = sumbech_optns_dct[submech_name]['fun_name'](
                species_list[0], self.spc_dct)
            self.species_subset_df = species_subset_df


        if len(species_list) > 0:
            self.mech_df_full = copy.deepcopy(self.mech_df)
            self.spc_dct_full = copy.deepcopy(self.spc_dct)

            self.mech_df, self.spc_dct = self.filter_byspecies(
                species_list, submech_name)
            self.species_list = species_list

    def filter_byspecies(self, species_list, filtertype):
        """ Find all reactions involving species of the species_list given as input
            Provides a new mechanism of the reactions of the selected subset
            and a reduced list of species

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.spc_dct: species dictionary
        :param species_list: list of species subsets
        :param filtertype: can be 'submech', 'submech_keepsubfuel/del', 'submech_prompt'

        :returns: mech_df, spc_dct
        :rtype: dataframe, dict
        """

        # add classification by subpes first
        if filtertype == 'submech_prompt':
            # reset the classification
            self.mech_df_full[['submech_prompt', 'rxn_ped']] = ''

        # deepcopy
        mech_df = copy.deepcopy(self.mech_df_full)
        # check that all species selected are in the species dictionary
        try:
            [self.spc_dct_full[sp] for sp in species_list]
        except KeyError as err:
            raise KeyError('Error in ISOLATE_SPECIES: ',
                  'not all species are in the species list. Exiting') from err

        # For all reactions in dataframe: check if species in selected list.
        # Otherwise remove the reaction
        spc_list = []
        for rxn in mech_df.index:
            rcts = list(mech_df['rct_names_lst'][rxn])
            prds = list(mech_df['prd_names_lst'][rxn])
            # Check if one species in the list is among reactants or products
            # of the reaction considered
            if filtertype in ['singlespecies', 'submech', 'submech_prompt']:
                _rchk = any(
                    rct == _spc for _spc in species_list for rct in rcts)
                _pchk = any(
                    prd == _spc for _spc in species_list for prd in prds)
                chk = int(_rchk or _pchk)

            elif filtertype in ['submech_deletelarge', 'submech_keepsubfuel']:

                _rchk = int(all(any(rct == _spc for _spc in species_list)
                               for rct in rcts))
                _pchk = int(all(any(prd == _spc for _spc in species_list)
                               for prd in prds))

                if filtertype == 'submech_deletelarge':
                    chk = int(_rchk and _pchk)
                    # all species of reactants and products have to be in the list
                elif filtertype == 'submech_keepsubfuel':
                    # filter out bimol/bimol reactions if some bimol rcts/prds are not in the species list
                    # equivalent to saying: at least all reactants (also single react works) or all products
                    # must be in the species list
                    chk = int(_rchk or _pchk)

            if chk >= 1 and filtertype != 'submech_prompt':
                spc_list.extend(rcts)
                spc_list.extend(prds)

            elif chk >= 1 and filtertype == 'submech_prompt':
                for sp in species_list:
                    if len(rcts) == 2 and any(rct == sp for rct in rcts) and len(prds) <= 2:
                        mech_df.at[rxn, 'submech_prompt'] = 'RAD_GEN_{}'.format(
                            sp)
                        mech_df.at[rxn, 'rxn_ped'] = '{}={}'.format(
                            rxn[0].split('=')[-1], rxn[0].split('=')[0])
                    elif len(prds) == 2 and any(prd == sp for prd in prds):
                        mech_df.at[rxn, 'submech_prompt'] = 'RAD_GEN_{}'.format(
                            sp)
                        mech_df.at[rxn, 'rxn_ped'] = rxn[0]

                    elif len(prds) > 2 and any(prd == sp for prd in prds):
                        mech_df.at[rxn, 'submech_prompt'] = 'PROMPT_LUMPED_{}'.format(
                            sp)
                    # rxn is unimol deco/formation of the radical
                    elif ((len(rcts) == 1 and rcts[0] == sp and len(prds) <= 2) or
                          (len(prds) == 1 and prds[0] == sp and len(rcts) <= 2)):
                        mech_df.at[rxn, 'submech_prompt'] = 'RAD_DECO_{}'.format(
                            sp)
                    elif (len(rcts) == 1 and rcts[0] == sp and len(prds) > 2):
                        mech_df.at[rxn, 'submech_prompt'] = 'RAD_DECO_LUMPED_{}'.format(
                            sp)

                    if len(mech_df.at[rxn, 'submech_prompt']) > 0:
                        break

            elif chk == 0 and filtertype != 'submech_prompt':  # reaction filtered out
                # don't filter for submech_prompt - you'll need it later to check for wellskipping channels
                mech_df = mech_df.drop(index=[rxn])

        if filtertype == 'submech_prompt':
            mech_df_new = pd.DataFrame(columns = mech_df.columns, dtype=object)
            # submech_prompt: if RAD_GEN/RAD_DECO in list, keep the rxns
            for _, subpes_df in mech_df.groupby(['pes', 'subpes']):
                if any('RAD' in CHECK or 'PROMPT' in CHECK   # if len > 1, value was assigned
                       for CHECK in subpes_df['submech_prompt'].values):
                    [spc_list.extend(list(rcts_tup))
                     for rcts_tup in subpes_df['rct_names_lst'].values]
                    [spc_list.extend(list(prds_tup))
                     for prds_tup in subpes_df['prd_names_lst'].values]
                    mech_df_new = pd.concat([mech_df_new, subpes_df], axis=0)
                    # add wellskipping channels that might be missing
                    if any('RAD_GEN' in CHECK
                           for CHECK in subpes_df['submech_prompt'].values):
                        added_rxns_df = self.add_wellskipping(subpes_df, mech_df_new)
                        mech_df_new = pd.concat([mech_df_new, added_rxns_df], axis=0)

            mech_df = copy.deepcopy(mech_df_new)

        # filter spc_list: unique elements
        spc_list = sorted(list(set(spc_list)))
        if filtertype == 'submech_deletelarge':
            musthaves = ['HE', 'AR', 'N2']
            [spc_list.append(m) for m in musthaves if m not in spc_list
            if m in self.spc_dct_full.keys()] # also consider 'must haves' that might not appear in reactions
        # new spc_dct
        spc_dct_val = list(map(self.spc_dct_full.get, spc_list))
        spc_dct = dict(zip(spc_list, spc_dct_val))

        return mech_df, spc_dct

    def add_wellskipping(self, subpes_df, mech_df_checkdup):
        """ check channels for a radical generation subset
            add temporary bimol or unimol wellskiping reactions generating the radical
            so they will be considered in the generation of prompt channels
            but removed in the final reaction list
            mech_df_checkdup is needed to check that the reaction is not for some reason already present in the mech
        """

        # identify radical(s) generating
        rad_list = []
        rad_bimol = []
        for rxn in sorted(subpes_df.index): # get sorted index
            if 'RAD_GEN' in subpes_df['submech_prompt'][rxn]:
                rad = subpes_df['submech_prompt'][rxn].split('_')[2]
                rad_list.append(rad)
                if rad in subpes_df['prd_names_lst_ord'][rxn]:
                    rad_bimol.append(subpes_df['prd_names_lst_ord'][rxn])
                elif rad in subpes_df['rct_names_lst_ord'][rxn]:
                    rad_bimol.append(subpes_df['rct_names_lst_ord'][rxn])

        rad_list = sorted(list(set(rad_list))) # reduce lists, might have found > 1 radical
        rad_bimol = sorted(list(set(rad_bimol)))

        # get reaction names
        rxn_list_ordered = list(
            zip(subpes_df['rct_names_lst_ord'].values, subpes_df['prd_names_lst_ord'].values))
        existing_rxn_list_ordered = list(
            zip(mech_df_checkdup['rct_names_lst_ord'].values, mech_df_checkdup['prd_names_lst_ord'].values))

        # get reactivity matrix for subpes
        connected_rxns_df = connect_rxn_df(rxn_list_ordered)
        # generate new reactions for the radicals
        new_wskip_rxns = []
        for rad_bim in rad_bimol:
            new_wskip_rxns.extend(add_wellskip(connected_rxns_df, rad_bim))

        new_wellskipping_idxs = []
        for rcts in new_wskip_rxns:
            # check that rxn is not already present in the inlet dataframe: if so, remove from list
            if tuple(rcts) in existing_rxn_list_ordered or tuple([rcts[1], rcts[0]]) in existing_rxn_list_ordered:
                # print('found {}'.format(rcts))
                continue
            # if not found, continue
            rxn = '{}={}'.format(
                '+'.join(rcts[0]), '+'.join(rcts[1]))
            new_wellskipping_idxs.append((rxn, (None,)))

        #print(subpes_df, '\n', rad_list, rad_bimol, '\n')
        wellskipp_rxns_df = pd.DataFrame(
            index=new_wellskipping_idxs, columns=subpes_df.columns, dtype=object)
        # common values
        wellskipp_rxns_df[['pes', 'subpes']
                          ] = subpes_df[['pes', 'subpes']].values[0]
        wellskipp_rxns_df['chnl'] = 'WELLSKIPPING FAKE'
        wellskipp_rxns_df['thrdbdy'] = [(None,)]*len(wellskipp_rxns_df.index)

        # add to dataframe
        for idx, rxn in enumerate(new_wellskipping_idxs):
            rad = sorted(list(set(rad_list).intersection(new_wskip_rxns[idx][1])))[0]

            wellskipp_rxns_df.at[rxn, 'submech_prompt'] = 'RAD_GEN_{}'.format(rad)
            wellskipp_rxns_df.at[rxn, 'rxn_ped'] = rxn[0]
            wellskipp_rxns_df.at[rxn, 'rct_names_lst'] = new_wskip_rxns[idx][0]
            wellskipp_rxns_df.at[rxn, 'prd_names_lst'] = new_wskip_rxns[idx][1]
            wellskipp_rxns_df.at[rxn, 'rct_names_lst_ord'] = new_wskip_rxns[idx][0]
            wellskipp_rxns_df.at[rxn, 'prd_names_lst_ord'] = new_wskip_rxns[idx][1]
            # print('ok for ', idx, rxn, rad_list, new_wskip_rxns[idx], '\n')

        return wellskipp_rxns_df

    def group_species(self, reac_sp_df):
        """ Checks if the reactions in self.mech_df contain any species
            of self.species_list and if so it marks the species.
            Assignment is hierarchical: the first species of species_list
            found determines the label of the reaction being checked

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.species_list: list of subset of species considered
        :param reac_sp_df: empty dataframe index=rxns, column: 'species'

        :returns: reac_sp_df dataframe[species][rxn]
        :rtype: dataframe[str][tuple]
        """

        # if species list is not found: do nothing, species entry remain empty
        if len(self.species_list) > 0:
            for rxn in reac_sp_df.index:
                rcts = list(self.mech_df['rct_names_lst'][rxn])
                prds = list(self.mech_df['prd_names_lst'][rxn])
                # check species hierarchically
                for sp_i in self.species_list:
                    if ((any(sp_i == rct for rct in rcts) or
                         any(sp_i == prd for prd in prds))
                            and isinstance(reac_sp_df['species'][rxn], float)):
                        reac_sp_df.at[rxn, 'species'] = sp_i
        return reac_sp_df

    def group_submech(self, submech_df):
        """ Assigns a submechanism (fuel, fuel radical, fuel add ..) to the reaction
            according to the species taking part to it.

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.species_list: list of subset of species considered
        :param self.species_subset_df: dataframe with species assigned to a
                                        certain type (fuel, fuel radical..)
        :param submech_df: empty dataframe index=rxns, column: 'submech'
                                                        or 'submech_keepsubfuel' or 'submech_deletelarge'
        :returns: submech_df dataframe with reactants
            mult dataframe[submech][rxn]
        :rtype: dataframe[str][tuple]
        """
        # better hierarchy:
        # if ANY of reactants belongs to
        # lbl as submech_df columns -> so it works with both submech and submech_keepsubfuel
        lbl_col = submech_df.columns[0]

        def assign_group(spcs_subset):
            # prioritize when a species is in fuel group;
            # assign "core" only when all species are in the core mech

            if any(subset in FUELGROUP for subset in spcs_subset):
                for fuelgroup in FUELGROUP:
                    if any(subset == fuelgroup for subset in spcs_subset):
                        group = fuelgroup
                        break
            elif all(subset in SUBMECHGROUP for subset in spcs_subset):
                for submechgroup in SUBMECHGROUP:
                    if any(subset == submechgroup for subset in spcs_subset):
                        group = submechgroup
                        break
            elif any(subset in SUPMECHGROUP for subset in spcs_subset):
                for supmechgroup in SUPMECHGROUP:
                    if any(subset == supmechgroup for subset in spcs_subset):
                        group = supmechgroup
                        break
            else:
                group = ''
                print('warning: species do not seem to belong to any group - check code')

            return group

        for rxn in submech_df.index:
            rcts = list(self.mech_df['rct_names_lst'][rxn])
            prds = list(self.mech_df['prd_names_lst'][rxn])
            spcs = rcts + prds
            # check if rcts and prds are in species list. if so, check spc type
            species_subset = []
            for spc in spcs:
                if spc in self.species_list:
                    species_subset.append(self.species_subset_df[spc])

            spcs_grp = assign_group(species_subset)
            # check species hierarchically (hierarchy fixed in species list)
            submech_df.at[rxn, lbl_col] = spcs_grp

        return submech_df

    def group_prompt(self, submech_df):
        """ group reactions according to prompt dissociation channels
            builds group dictionary list:
            :param grps: [{'grp': 0, 'idxs': [], 'peds': [], 'hot': []}, ...]
            :type grps: list(dct)
        """

        # assign group values if not present:
        # 1. groups by PES
        # 2. if PES is rad_deco: renames as rad_deco_radname
        # 3. if PES produces hot radical: adds it to a group ['grp N']
        #       then specifies the name in ['groupname'] nb add to same list if same subpes
        try:
            self.grps
        except AttributeError:
            grps = []
            self.species_deco_dct = dict.fromkeys(self.species_list)
            # filter by RAD_DECO: save the pes/subpes value

            for hot_label in ['RAD_DECO_{}'.format(sp) for sp in self.species_list]:
                hot_df = self.mech_df[self.mech_df['submech_prompt']
                                      == hot_label]
                sp = hot_label.split('_')[-1]
                self.species_deco_dct[sp] = '{}:{}'.format(
                    hot_df['pes'].iloc[0], hot_df['subpes'].iloc[0])

            grp_dct_template = {'grp': 0, 'idxs': [], 'peds': [], 'hot': []}
            grpN = 0
            rad_gen_list = ['RAD_GEN_{}'.format(
                sp) for sp in self.species_list]
            ped_df = self.mech_df[self.mech_df['submech_prompt'].isin(
                rad_gen_list)]

            for pes, pesdf in ped_df.groupby('pes'):
                grpN += 1
                hotspc_dct = {}  # 'pes:subpes':[sp1, sp2]
                grp_dct = copy.deepcopy(grp_dct_template)
                grp_dct['grp'] = grpN

                for subpes, subpesdf in pesdf.groupby('subpes'):
                    peds = []
                    grp_dct['hot'].append([])
                    grp_dct['idxs'].append('{}:{}'.format(pes, subpes))
                    rxns = subpesdf.sort_values(by=['rxn_names'], ascending=False).index #ordered like the channels
                    for rxn in rxns:
                        peds.append(subpesdf['rxn_ped'][rxn])
                        # deal with hotspecies
                        sp = subpesdf['submech_prompt'][rxn].split('_')[-1]
                        hotpes = self.species_deco_dct[sp]
                        if hotpes in hotspc_dct.keys():
                            hotspc_dct[hotpes].append(sp)
                        else:
                            hotspc_dct[hotpes] = [sp]
                    grp_dct['peds'].append(peds)
                # assign hot species
                for hotpes in hotspc_dct.keys():
                    grp_dct['peds'].append([])
                    grp_dct['idxs'].append(hotpes)
                    grp_dct['hot'].append(sorted(list(set(hotspc_dct[hotpes]))))
                    grp_dct['modeltype'] = 'rovib_dos'

                grps.append(grp_dct)

            self.grps = grps

        submech_df = pd.DataFrame(
            index=self.mech_df.index)

        return submech_df # aligned with other functions, but does not return a dataframe

    def filter_groups_prompt(self, therm_dct, DFG, T0=300.):
        """ filter prompt groups according to DFG dct criteria:
            A+B->B+C->B+(D+E) must be < X kcal/mol endothermic
            If DH1 is added to dissociation reaction: k(T*)/k(T) > RATIO
            ratio: (DH1+DH3)/DH3 < threshold
            Absolute threshold for k(T*): should be > k*abs (only for exothermic rxns not fulfilling other criteria)
            NB all criteria should be strict --> by default, no reaction is kept.
            the criterion introduced should be the strict one
            T0 is the T at which DH are considered
            Tref is the T at which the decomposition rate of the radical is extracted
            otherwise remove from dct
        """

        # define criteria
        self.DHmax = DFG['DH']
        self.H5H3ratio = DFG['H5H3ratio']
        self.kratio = DFG['kratio']
        self.kabs = DFG['kabs']
        self.Tref = DFG['Tref']
        self.keepfiltered = int(DFG['keepfiltered'])
        self.lookforpromptchains = bool(DFG['lookforpromptchains'])
        # convert therm dct to dataframe
        self.therm_df = thermo.spc_therm_dct_df(therm_dct)
        # 1 find min endothermicity of deco reactions of hotspecies
        self.dh_min_hot = dict.fromkeys(self.species_list)
        self.k_max_hot = dict.fromkeys(self.species_list)
        self.labels_hot = dict.fromkeys(self.species_list)
        hot_sp_df_dct = dict.fromkeys(self.species_list)


        for hot_sp in self.species_list:
            hot_sp_df_dct[hot_sp] = self.mech_df[self.mech_df['submech_prompt']
                                                 == 'RAD_DECO_{}'.format(hot_sp)]
            self.k_max_hot[hot_sp], self.dh_min_hot[hot_sp], self.labels_hot[hot_sp] = nonboltz.get_max_reactivity(
                hot_sp, hot_sp_df_dct[hot_sp], self.therm_df, T0, self.Tref)

        # loop over groups and delete too endothermic reactions
        filtered_grps = []
        grp_idx = 0
        fmt_lbls = numpy.array(['%40s', '%15s', '%.1f', '%.1f',
                                '%.1f', '%1.1e', '%1.1e', '%.0f', '%s'], dtype=object)
        lbls = numpy.array(['prompt rxn{}\t'.format(' '*40), 'hot species', 'dh1({:.0f}K)*phi'.format(T0), 'dh2({:.0f}K)'.format(T0),
                            'dhtot({:.0f}K)'.format(T0), 'k({:.0f}K)'.format(self.Tref), 'k*(T*({:.0f}K))'.format(self.Tref), 'T*({:.0f}K)'.format(self.Tref), 'keep?'], dtype=object)

        self.rxns_dh = numpy.vstack((fmt_lbls, lbls))


        for grp in self.grps:
            grp_new = {'grp': 0, 'idxs': [],
                       'peds': [], 'hot': [], 'modeltype': ''}
            self.grp_add = {'idxs': [], 'peds': [], 'hot': []}
            check = 0
            check_filtered = 0
            active_hotsp = []
            # lists, potentially more than 1
            for n, ped in enumerate(grp['peds']):
                exceptions = 0
                newped = []
                for ped_i in ped:
                    check_ped_i = 0
                    # extract rxn exo/endo thermicity
                    print('analyze ped .. {}'.format(ped_i))
                    rcts = ped_i.split('=')[0].split('+')
                    prds = ped_i.split('=')[1].split('+')

                    hot_spcs = sorted(list(set(self.species_list) & set(prds)))

                    for hot_sp in hot_spcs:

                        nonhot = list(set([hot_sp]) ^ set(prds))
                        if len(nonhot) == 0:
                            nonhot = hot_sp #it means you have something like A+B=>2C (disproport. or dissociation)
                        else:
                            nonhot = nonhot[0]

                        phi = phi_equip_fromdct(hot_sp, nonhot, self.spc_dct)
                        try:
                            dh = thermo.extract_deltaX_therm(
                                self.therm_df, rcts, prds, 'H')/1000
                        except KeyError as keyerr:
                            # it is possible that the species has no thermo, e.g., unstable species.
                            print('thermo not found for:', keyerr.args[0])
                            exceptions += 1
                            # save hotsp anyway since you don't know what's going to happen
                            active_hotsp.append(hot_sp)
                            continue
                        # calculate equivalent T* at Tref and corresponding k* rate
                        # phi: fraction of en transferred to prods - very approixmate
                        try:
                            T_star, k_star, dh_tot = nonboltz.estimate_hot_hk(
                                dh*phi*1000, self.Tref, self.therm_df[hot_sp]['Cp'], self.k_max_hot[hot_sp], self.dh_min_hot[hot_sp]*1000)
                        except TypeError:
                            print('dh failed for:', rcts, prds, hot_sp)
                            exceptions += 1
                            # save hotsp anyway since you don't know what's going to happen
                            active_hotsp.append(hot_sp)
                            continue

                        # print(dh_tot[T0], dh_tot[T0]/self.dh_min_hot[hot_sp][T0]*int(self.dh_min_hot[hot_sp][T0]), k_star/self.k_max_hot[hot_sp][self.Tref], k_star, dh[T0], self.dh_min_hot[hot_sp][T0])
                        # print(dh_tot[T0] < self.DHmax, dh_tot[T0]/self.dh_min_hot[hot_sp][T0]*int(self.dh_min_hot[hot_sp][T0] > 0) < self.H5H3ratio, k_star/self.k_max_hot[hot_sp][self.Tref] > self.kratio, k_star > self.kabs, dh[T0] < 0)

                        if dh_tot[T0] < self.DHmax or dh_tot[T0]/self.dh_min_hot[hot_sp][T0]*int(self.dh_min_hot[hot_sp][T0] > 0) < self.H5H3ratio \
                            or (k_star/self.k_max_hot[hot_sp][self.Tref] > self.kratio) \
                                or (k_star > self.kabs and dh[T0] < 0):
                            check += 1
                            keep = 'YES'
                            active_hotsp.append(hot_sp)
                            # BUILD REACTION CHAINS STARTING FROM HERE
                            # look for chains
                            if self.lookforpromptchains == True:
                                self.rxn_chain_prompt(
                                    T0, dh_tot, hot_sp, hot_sp_df_dct)
                            check_ped_i += 1
                        else:
                            if self.keepfiltered != 0:
                                check_ped_i += 1
                                active_hotsp.append(hot_sp)
                                check += 1
                                check_filtered += 1  # count filtered rxns
                            keep = 'NO'

                        array_info = numpy.array(
                            [ped_i, hot_sp, dh[T0]*phi, self.dh_min_hot[hot_sp][T0], dh_tot[T0], self.k_max_hot[hot_sp][self.Tref], k_star, T_star, keep], dtype=object)
                        self.rxns_dh = numpy.vstack((self.rxns_dh, array_info))

                    if check_ped_i >= 1:
                        newped.append(ped_i)

                if exceptions == len(newped) and len(newped) > 0:
                    check = 1  # keep things you were unable to compute stuff for that ped
                    print(
                        'Warning: unable to derive thermo / rates for set of peds: {}'.format(newped))
                    print('these peds will be kept')

                if newped:
                    grp_new['idxs'].append(grp['idxs'][n])
                    grp_new['peds'].append(newped)
                    grp_new['hot'].append([])
                    # CHECK FOR NEW CHAINS AND UPDATE GRP NEW
                    # NB REMOVE hotsp FROM ACTIVE_HOTSP TO AVOID DOUBLE COUNTING

            for n, hot in enumerate(grp['hot']):
                grp['hot'][n] = [hot_i for hot_i in hot if hot_i in active_hotsp]
                grp['hot'][n].sort()  # keep same order
                if grp['hot'][n]:
                    grp_new['idxs'].append(grp['idxs'][n])
                    # APPEND ALSO HOT SPECIES IF NEW CHAINS CONSIDERED
                    grp_new['peds'].append([])
                    grp_new['hot'].append(grp['hot'][n])

            # these are single indices
            for n, idx in enumerate(self.grp_add['idxs']):
                # merge other groups -
                # if you have same idxs: it is ped + hot, so take the
                # hot of the grp_new (more complete) and the ped label of the latter
                try:
                    i_idx = grp_new['idxs'].index(idx)
                    grp_new['peds'][i_idx] = self.grp_add['peds'][n]
                except ValueError:
                    grp_new['idxs'].append(idx)
                    # APPEND ALSO HOT SPECIES IF NEW CHAINS CONSIDERED
                    grp_new['peds'].append(self.grp_add['peds'][n])
                    grp_new['hot'].append(self.grp_add['hot'][n])

            if check > 0:
                grp_idx += 1
                grp_new['grp'] = grp_idx
                if check_filtered == check:
                    grp_new['modeltype'] = 'thermal'
                else:
                    grp_new['modeltype'] = 'rovib_dos'
                filtered_grps.append(grp_new)

        self.grps = filtered_grps

        # remove fake rxns - unneeded because you already do in sort_and_label
        # rxns_fake = self.mech_df[self.mech_df['chnl']
        #                          == 'WELLSKIPPING FAKE'].index
        # self.mech_df = self.mech_df.drop(index=rxns_fake)
        # resort because you added reactions
        self.sort_and_label(self.criteria_all, self.labels_all)

    def rxn_chain_prompt(self, T0, dh_tot, rad, sp_df_dct):
        """ from a given hot product (rad), derive prompt reaction chain complying with thresholds
        """
        # append to self groups only if pes not present as ped
        if self.species_deco_dct[rad] not in self.grp_add['idxs']:
            self.grp_add['idxs'].append(self.species_deco_dct[rad])
            self.grp_add['peds'].append([])
            self.grp_add['hot'].append([rad])

        check_break = 0
        while check_break == 0:
            # radical generating reaction (and hot sp)

            dh_start = copy.deepcopy(dh_tot)
            # rad is the potentially prompt species:
            # 1. CHECK decomposition channels:
            # 1.1 DECO CHANNEL: IF IT PRODUCES RADICALS, CONTINUE
            prds = self.labels_hot[rad].split('=')[1].split('+')
            # potentially change the algorithm:
            # if you have 2 products, just analyze both; select the one
            # most 'compliant' with the conditions (e.g., smallest dhtot ratio)
            # and go on with that.

            if len(prds) != 2:
                break

            phis, T_stars, k_stars, dh_tots, add_check, hot_mech_df, hot_spc_dct, pesN, subpesN = (
                {prd: None for prd in prds} for _ in range(9))

            for prd in prds:
                if sum(automol.chi.formula(
                        self.spc_dct[prd]['inchi']).values()) < 3:
                    continue

                if prd not in self.k_max_hot.keys():
                    add_check[prd] = 1
                    hot_mech_df[prd], hot_spc_dct[prd] = self.filter_byspecies(
                        [prd], 'submech_prompt')
                    sp_df_dct[prd] = hot_mech_df[prd][hot_mech_df[prd]['submech_prompt']
                                                 == 'RAD_DECO_{}'.format(prd)]

                    if len(sp_df_dct[prd]) == 0:
                        continue  # it's possible that the radical deco is not present in the mech!
                    pesN[prd] = sp_df_dct[prd]['pes'].iloc[0]
                    subpesN[prd] = sp_df_dct[prd]['subpes'].iloc[0]
                    self.species_deco_dct[prd] = '{}:{}'.format(pesN[prd], subpesN[prd])
                    self.k_max_hot[prd], self.dh_min_hot[prd], self.labels_hot[prd] = nonboltz.get_max_reactivity(
                        prd, sp_df_dct[prd], self.therm_df, T0, self.Tref)  # high T to get bimol faster
                else:
                    add_check[prd] = 0

                try:
                    nonprd = list(set([prd]) ^ set(prds))[0]
                    phis[prd] = phi_equip_fromdct(prd, nonprd, self.spc_dct)

                    T_stars[prd], k_stars[prd], dh_tots[prd] = nonboltz.estimate_hot_hk(
                        dh_start*phis[prd]*1000, self.Tref, self.therm_df[prd]['Cp'], self.k_max_hot[prd], self.dh_min_hot[prd]*1000)

                except TypeError:
                    print('dh failed for: {}'.format(self.labels_hot[rad]))
                    break

            # identify hot species by min dh_tot
            if not any(phis.values()):
                break
            elif not all(phis.values()):
                hot_sp = [key for key in phis.keys() if phis[key]][0]
            else:
                dh_tots_T0 = [dh_tots[prds[0]][T0], dh_tots[prds[1]][T0]]
                min_dhtot = min(dh_tots_T0)
                hot_sp = prds[dh_tots_T0.index(min_dhtot)]

            dh_tot = dh_tots[hot_sp]
            print('prompt chain for {}'.format(self.labels_hot[rad]),
                  'hot species is {}'.format(hot_sp),
                  'fraction of en transferred: {:.1f} \n'.format(phis[hot_sp]))

            ##################################################
            if add_check[hot_sp] == 1:
                # # UPDATE SELF.MECH_DF AND SPECIES TO ADD NEW PESs
                # update mechanism: include all elements in the subpes of interest
                # print(hot_sp, hot_mech_df[hot_sp])
                new_hot_mech_df = hot_mech_df[hot_sp][(hot_mech_df[hot_sp]['pes'] == pesN[hot_sp]) & (hot_mech_df[hot_sp]['subpes'] == subpesN[hot_sp])]
                self.mech_df = pd.concat(
                    [self.mech_df, new_hot_mech_df], axis=0)
                # print(new_hot_mech_df)
                # update species avoiding duplicate entries
                new_spc = [sp for sp in hot_spc_dct[hot_sp].keys() if sp not in self.spc_dct.keys()]
                self.spc_dct.update(dict(zip(new_spc, list(map(hot_spc_dct[hot_sp].get, new_spc)))))

            # 2.2 IF COMPLIANT WITH CONDITIONS:
            if dh_tot[T0] < self.DHmax or dh_tot[T0]/self.dh_min_hot[hot_sp][T0]*int(self.dh_min_hot[hot_sp][T0] > 0) < self.H5H3ratio \
                or (k_stars[hot_sp]/self.k_max_hot[hot_sp][self.Tref] > self.kratio) \
                    or (k_stars[hot_sp] > self.kabs and dh_start[T0] < 0):
                keep = 'YES'
                array_info = numpy.array(
                    [self.labels_hot[rad], hot_sp, dh_start[T0]*phis[hot_sp], self.dh_min_hot[hot_sp][T0],
                     dh_tot[T0], self.k_max_hot[hot_sp][self.Tref], k_stars[hot_sp], T_stars[hot_sp], keep], dtype=object)
                self.rxns_dh = numpy.vstack((self.rxns_dh, array_info))
                peds = [self.labels_hot[rad]]

            else:
                # keep species only as hot
                check_break = 1
                peds = []

            if peds != []:
                # find index and replace ped
                idx = self.grp_add['idxs'].index(self.species_deco_dct[rad])
                self.grp_add['peds'][idx] = peds

                if self.species_deco_dct[hot_sp] not in self.grp_add['idxs']:
                    self.grp_add['idxs'].append(self.species_deco_dct[hot_sp])
                    self.grp_add['peds'].append([])
                    self.grp_add['hot'].append([hot_sp])

            # updates for next cycle
            #     - new cycle: rad = hot_sp

            rad = copy.deepcopy(hot_sp)

    def reac_mult(self, reac_mult_df):
        """ determines the multiplicity of the reactants of the reactions
            considered

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.spc_dct: species dictionary
        :param rxncl_mult_df: empty dataframe index=rxns, column: 'mult'

        :returns: reac_mult_df dataframe with reactants mult
            dataframe[mult][rxn]
        :rtype: dataframe[str][tuple]
        """
        # assign multiplicity values to each reactant
        for rxn in reac_mult_df.index:
            mult = 1
            rcts = self.mech_df['rct_names_lst'][rxn]
            mult = get_mult(rcts, self.spc_dct)
            reac_mult_df.at[rxn, 'mult'] = str(mult)

        return reac_mult_df

    def rxnclass_broad(self, rxncl_broad_df):
        """ assigns reaction classes based on stoichiometry and molecularity

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param self.spc_dct: species dictionary
        :param rxncl_broad_df: empty dataframe
            index=rxns, column: 'rxn_class_broad'

        :returns: rxncl_broad_df dataframe with reaction classes
            dataframe[class][rxn]
        :rtype: dataframe[str][tuple]
        """
        for rxn in rxncl_broad_df.index:
            rcts = self.mech_df['rct_names_lst_ord'][rxn]
            prds = self.mech_df['prd_names_lst_ord'][rxn]
            _smol = (self.mech_df['molecularity'][rxn] == 1)
            _tbody = (self.mech_df['molecularity'][rxn] == 2
                      and self.mech_df['isthrdbdy'][rxn] == 1)

            if _smol or _tbody:
                # unimolecular reaction classification
                rxn_class_broad = rxnclass.classify_unimol(
                    rcts, prds, self.spc_dct)
            else:
                # bimolecular reaction classification
                rxn_class_broad = rxnclass.classify_bimol(
                    rcts, prds, self.spc_dct)
            rxncl_broad_df.at[rxn, 'rxn_class_broad'] = rxn_class_broad
            # (rxncl_broad_df.loc[[rxn], 'rxn_class_broad']) # technicall

        return rxncl_broad_df

    def rxnclass_graph(self, rxncl_graph_df):
        """ assigns reaction classes using graph approach to all reactions
            first subdivides the mech into subpeses; then classifies all rxn
            within the subpes, including wellskipping channels

        :param self.mech_df: dataframe with mech info (contains all reactions)
        :param rxncl_graph_df: empty dataframe
            index=rxns, column: 'rxn_class_graph'

        :returns: rxncl_graph_df dataframe with reaction classes
            dataframe[class][rxn]
        :rtype: dataframe[str][tuple]
        """

        # 1. Group by subpes if absent
        # done in __init__ in the updated version
        # self.mech_df = pd.concat([self.mech_df, self.chnl('')], axis=1)

        # 2. Graph classification or each subpes
        for _, subpes_df in self.mech_df.groupby(['pes', 'subpes']):
            # sort by molecularity: analyze first unimolecular isomerizations,
            # unimolecular decompositions, and then bimolecular reactions
            subpes_df = subpes_df.sort_values(
                by=['molecularity', 'N_of_prods'])
            # REFER TO REORDERED SPECIES NAMES,
            # OTHERWISE YOU MAY HAVE INCONSISTENT SPECIES NAMING
            # subpes species list
            species_subpes = sorted(list(set(
                list(subpes_df['rct_names_lst_ord'].values) +
                list(subpes_df['prd_names_lst_ord'].values))))

            # build dataframe: elementary reactivity matrix
            _nsubpes = len(species_subpes)
            elem_reac_df = pd.DataFrame(
                numpy.zeros((_nsubpes, _nsubpes), dtype='<U32'),
                index=species_subpes,
                columns=species_subpes)

            # graph classification
            for rxn in subpes_df.index:
                rct_names = subpes_df['rct_names_lst'][rxn]
                prd_names = subpes_df['prd_names_lst'][rxn]
                rct_names_ord = subpes_df.at[rxn, 'rct_names_lst_ord']
                prd_names_ord = subpes_df.at[rxn, 'prd_names_lst_ord']
                # Exclude rxns with more than 2 rcts or prds (not elementary!)
                if len(rct_names) < 3 and len(prd_names) < 3:

                    rclass = rxnclass.classify_graph(
                        self.spc_dct, rct_names, prd_names)

                else:
                    rclass = 'unclassified - lumped'
                rxncl_graph_df.at[rxn, 'rxn_class_graph'] = rclass

                # store values in the elementary reactivity matrix
                # (for now contaminated with isomerizations)
                elem_reac_df.at[prd_names_ord, rct_names_ord] = rclass
                elem_reac_df.at[rct_names_ord, prd_names_ord] = rclass

            # 3. classify well skipping channels
            # reclassify the unclassified reactions A->B+C, B+C->D, B+C->E+F
            for rxn in subpes_df.index:
                if rxncl_graph_df['rxn_class_graph'][rxn] == 'unclassified':

                    # call external function for WS channel classification

                    rxn_type_ws = rxnclass.classify_ws(
                        subpes_df, elem_reac_df, species_subpes, rxn)
                    if rxn_type_ws is not None:
                        rxncl_graph_df.at[rxn, 'rxn_class_graph'] = rxn_type_ws

        return rxncl_graph_df

    def rxn_max_vals(self, rxn_maxvals_df):
        """ Determines the maximum value of the rates in ktp dictionary.

        :param rxn_maxvals_df:
            empty dataframe index=rxns, column: 'rxn_max_vals'

        :returns: rxn_maxvals_df dataframe with overall max value of the rate
        :rtype: dataframe[float][tuple]
        """
        # extract maximum value for each ktp dictionary
        for rxn in rxn_maxvals_df.index:
            param_vals_dct = self.mech_df['param_vals'][rxn]
            max_val = calc_rates.get_max_aligned_values(param_vals_dct)
            rxn_maxvals_df.at[rxn, 'rxn_max_vals'] = max_val

        return rxn_maxvals_df

    def rxn_max_ratio(self, rxn_maxratio_df):
        """ Determines the maximum value of the ratios between
            different rates of ktp dct.
            (compares two dictionaries at the same pressure)

        :param rxn_maxratio_df:
            empty dataframe index=rxns, column: 'rxn_max_ratio'

        :returns: rxn_max_ratio dataframe with overall max value of the ratio
        :rtype: dataframe[float][tuple]
        """
        # extract maximum ratio for each set ktp dictionary
        for rxn in rxn_maxratio_df.index:
            param_vals_dct = self.mech_df['param_vals'][rxn]
            # get the ratio:
            param_ratio_dct = calc_rates.get_aligned_rxn_ratio_dct(
                param_vals_dct)
            max_val = calc_rates.get_max_aligned_values([param_ratio_dct[1]])
            rxn_maxratio_df.at[rxn, 'rxn_max_ratio'] = max_val

        return rxn_maxratio_df

    # ASSIGN HEADERS #

    def class_headers(self, hierarchy, labels):
        """ assigns class headers based on the hierarchy provided in the input


        :param self.mech_df: dataframe with mech info
        :param hierarchy: list of strings with hierarchical order
            of sorting criteria ['..','..','N']: N is of criteria for headers
            all the other criteria are written as inline comments
        :param labels: labels corresponding to all possible sorting criteria
                        labels are then written as comments

        :returns: None. updates self.mech_df['cmts_top','cmts_inline']
                    'cmts_top': comments to write as header for a reaction
                    'cmts_inline': comments to write on same line of reaction
        """
        # if comments already in mech_df: remove them:
        if 'cmts_inline' in self.mech_df.columns and 'cmts_top' in self.mech_df.columns:
            self.mech_df = self.mech_df.drop(
                ['cmts_inline', 'cmts_top'], axis=1)
        # df for comments_top and comments_inline
        ept_df = numpy.zeros((len(self.mech_df.index), 1), dtype=str)
        df_cmts_top = pd.DataFrame(
            ept_df, index=self.mech_df.index, columns=['cmts_top'])
        df_cmts_inline = pd.DataFrame(
            ept_df, index=self.mech_df.index, columns=['cmts_inline'])

        try:
            n_headers = int(hierarchy[-1])
        except ValueError as err:
            raise ValueError(
                '*ERROR: Last line of sorting options is the N ',
                'of criteria to use for class headers') from err

        # Write topheader comments
        if n_headers > 0:
            for name, rdf in self.mech_df.groupby(hierarchy[:n_headers]):
                # Write rxn class as top header comments
                rxnclass = cmts_string(
                    name, labels[hierarchy[:n_headers]], 'class_head')
                idx0 = rdf.index[0]
                df_cmts_top.at[idx0, 'cmts_top'] = rxnclass

                # Write inline comments if necessary
                if n_headers < len(hierarchy)-1:
                    for name2, rdf2 in rdf.groupby(hierarchy[n_headers:-1]):
                        rxnclass = cmts_string(
                            name2, labels[hierarchy[n_headers:-1]], 'subclass')
                        df_cmts_inline.loc[rdf2.index, 'cmts_inline'] = rxnclass
        else:
            # Write only inline comments
            for name, rdf in self.mech_df.groupby(hierarchy[n_headers:-1]):
                rxnclass = cmts_string(
                    name, labels[hierarchy[n_headers:-1]], 'class')
                df_cmts_inline.loc[rdf.index, 'cmts_inline'] = rxnclass

        # Concatenate DFs
        self.mech_df = pd.concat(
            [self.mech_df, df_cmts_top, df_cmts_inline], axis=1)

    # OUTPUT DATAFRAMES and DICTIONARIES #
    def return_mech_df(self):
        """ Provides sorted rxn indices and associated comments, as well as
            sorted species dictionary.

        :param self.mech_df: dataframe with mech info
        :param self.spc_dct: species dictionary for reactions of sorted mech

        :returns: sorted reaction names, comments, sorted species dictionary
        :rtype: list[tuple],
                dict{tuple(r1, r2,): {cmts_top: str, cmts_inline: str}, ...},
                dict
        """

        rct_names = self.mech_df['rct_names_lst'].values
        prd_names = self.mech_df['prd_names_lst'].values
        thrdbdy = self.mech_df['thrdbdy'].values
        new_idx = list(zip(rct_names, prd_names, thrdbdy))
        # store comments in dct
        cmts_df = pd.DataFrame(
            self.mech_df[['cmts_top', 'cmts_inline']].values,
            index=new_idx,
            columns=['cmts_top', 'cmts_inline'])
        cmts = cmts_df.to_dict('index')

        return new_idx, cmts, self.spc_dct

    def return_pes_dct(self):
        """ returns a PES dictionary

            chn idx is set to the OVERALL PES, not the SUB_PES

        :returns: pes_dct: {(fml, n_pes, n_subpes): ((chn_idx, (rcts, prds)),)}
        :rtype: dct{tuple: tuple}
        """

        # get the pes dictionary
        pes_dct = {}

        for pes_idx, pes_all_df in self.mech_df.groupby('pes'):
            fml_str = pes_all_df['formulas'].values[0]

            for subpes_idx, pes_dct_df in pes_all_df.groupby('subpes'):
                # Get the ('fml', n_pes, n_subpes) for dict key
                pes_dct_df = pes_dct_df.sort_values(by='chnl')
                pes_dct_key = (fml_str, pes_idx-1, subpes_idx-1)
                pes_chnl_tuples = pes_dct_df['pes_chnl_tuple'].values
                chnl_lst = tuple(pes_chnl_tuples)
                pes_dct[pes_dct_key] = chnl_lst

        return pes_dct
