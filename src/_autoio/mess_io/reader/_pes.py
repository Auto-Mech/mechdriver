"""
  Read a MESS input file and compile data for the PES inside
"""

import numpy as np
import autoparse.find as apf
import autoparse.pattern as app
from ioformat import remove_comment_lines


def pes(input_string, read_fake=False):
    """ Read a MESS input file string and get info about PES

        :param input_string: string for a MESS (rates) input file
        :type input_string: str
        :param read_fake: value to include fake wells and barriers
        :type read_fake: bool
        :return energy_dct: dict[label: energy]
        :rtype: dict[label: energy]
        :return conn_lst_dct: defines the wells connected by each barrier
        :rtype: dict[barrier: (reac, prod)] all str
        :return conn_lst
        :rtype: lst(str)
    """

    # Initialize energy and connection information
    energy_dct = {}
    conn_lst = tuple()
    conn_lst_dct = {}
    pes_label_dct = {}

    input_string = remove_comment_lines(
        input_string, delim_pattern=app.escape('!'))
    input_lines = input_string.splitlines()
    for idx, line in enumerate(input_lines):

        if 'Well' in line:

            line_lst = line.strip().split()
            if line_lst[0].strip() == 'Well' and '!' not in line_lst[0]:
                # Get label
                label = line_lst[1]

                if ('F' not in label) or ('F' in label and read_fake):
                    # Get energy
                    for line2 in input_lines[idx:]:
                        if 'ZeroEnergy' in line2:
                            ene = float(line2.split()[-1])
                            break

                    # Add value to energy dct
                    energy_dct[label] = ene

                    line_lst2 = line.split('!')
                    if len(line_lst2) == 1:
                        line_lst2 = line.split('#')
                    if len(line_lst2) > 1:
                        spc = line_lst2[1]
                        pes_label_dct[spc.strip()] = label
                    else:
                        pes_label_dct[label] = label
                        #print('Warning: labeling not found for '
                        #      f'species {label}')

        if 'Bimolecular' in line:

            line_lst = line.strip().split()

            if line_lst[0].strip() == 'Bimolecular' and '!' not in line:
                # Get label
                label = line_lst[1]

                # Get energy
                for line2 in input_lines[idx:]:
                    if 'Dummy' in line2:
                        ene = -10.0
                        break
                    if 'GroundEnergy' in line2:
                        ene = float(line2.split()[-1])
                        break

                # Add value to dct
                energy_dct[label] = ene

                # Add value to PES dct - NB THIS DEPENDS ON THE INPUT FILE.
                # IF NOT PRESENT, DO NOT GENERATE THE PES LABEL DICTIONARY
                cnt = 0
                frags = []
                for line2 in input_lines[idx:]:
                    if 'Fragment' in line2:
                        # Try and grab name from comment line
                        frag_line_lst = line2.split('!')
                        if len(frag_line_lst) == 1:
                            frag_line_lst = line2.split('#')
                        if len(frag_line_lst) > 1:
                            frag = frag_line_lst[1]
                            # strip gets rid of the spaces before and after
                            frags.append(frag.strip())
                        else:
                            frag_line_lst = line2.split()
                            frag = frag_line_lst[1]
                            frags.append(frag.strip())
                            #print('Warning: labeling not found for '
                            #      f'bimol fragments for {label}')
                        cnt += 1
                    if cnt == 2:
                        break

                pes_label_dct[' + '.join(frags)] = label

        if 'Barrier' in line:

            line_lst = line.strip().split()
            if line_lst[0].strip() == 'Barrier' and '!' not in line:
                # Get label
                [tslabel, rlabel, plabel] = line_lst[1:4]

                if ('F' not in tslabel) or ('F' in tslabel and read_fake):
                    # Get energy
                    for line2 in input_lines[idx:]:
                        if 'ZeroEnergy' in line2:
                            ene = float(line2.split()[-1])
                            break

                    # Add value to dct
                    energy_dct[tslabel] = ene

                    # Amend fake labels (may be wrong)
                    if not read_fake:
                        rlabel = rlabel.replace('F', 'P')
                        plabel = plabel.replace('F', 'P')

                    # Add the connection to lst
                    conn_lst += ((rlabel, tslabel),)
                    conn_lst += ((tslabel, plabel),)
                    conn_lst_dct[tslabel] = (rlabel, plabel)

    return energy_dct, conn_lst, conn_lst_dct, pes_label_dct


def find_barrier(conn_lst_dct, reac, prod):
    """ finds the barrier that connects reac to prod
        returns None if the barrier is not found
        future implementation: should find lowest energy path from reac to prod

        :param conn_lst_dct: defines the wells connected by each barrier
        :type conn_lst_dct: dict[barrier: (reac, prod)] all str
        :param reac, prod: connected species
        :type reac, prod: str
        :return barriername: name of the barrier
        :rtype: str
    """

    bar1 = (reac, prod)
    bar2 = (prod, reac)
    ls_keys = list(conn_lst_dct.keys())
    ls_vals = list(conn_lst_dct.values())
    bar1_find = apf.where_is(bar1, ls_vals)
    bar2_find = apf.where_is(bar2, ls_vals)
    if len(bar1_find) > 0:
        _bar = ls_keys[bar1_find[0]]
    elif len(bar2_find) > 0:
        _bar = ls_keys[bar2_find[0]]
    else:
        _bar = None

    return _bar


def get_species(input_string):
    """ Read a MESS input file string and get the block of each species
        Bimolecular fragments are listed together,
        but header Fragment is changed to Species

        :param input_string: string for a MESS (rates) input file
        :type input_string: str

        :return species_blocks: dictionary with the species blocks
            {name:[frag1 block, frag2 block], name:[unimol block],}
        :rtype: dict{label: list}
    """
    # find where data of interest are
    bad_wrds = ['WellDepth', 'WellCutoff', 'WellExtension',
                    'WellReductionThreshold', 'WellPartitionMethod',
                    'WellProjectionThreshold', 'PEDSpecies']
    
    input_string = remove_comment_lines(
        input_string, delim_pattern=app.escape('!'))
    lines = input_string.splitlines()
    lines = [line for line in lines if (line.strip() != '' and
             all(bad not in line for bad in bad_wrds))]

    _name_arr = np.array(
        [('Bimolecular' in line or 'Well' in line)
         for line in lines],
        dtype=int)
    names_i = np.where(_name_arr == 1)[0]

    _init_arr = np.array(
        [('Fragment' in line or 'Species' in line) and
         'FragmentGeometry' not in line for line in lines],
        dtype=int)
    init_i = np.where(_init_arr == 1)[0]
    init_i = init_i[init_i > names_i[0]]
    
    _barriers = np.array(['Barrier' in line for line in lines], dtype=int)
    first_barrier = np.where(_barriers == 1)[0]
    final_i = np.append(init_i[1:], first_barrier) -1 # line before each name

    # dictionary labels
    labels = [lines[i].strip().split()[1] for i in names_i]
    species_blocks = {k: [] for k in labels}

    # extract the data
    for i in np.arange(0, len(init_i)):

        # type
        sp_type = lines[init_i[i]].strip().split()[0]
        label = lines[names_i[init_i[i] > names_i][-1]].strip().split()[1]

        # name and label
        if sp_type == 'Fragment':
            name = 'Species ' + lines[init_i[i]].strip().split()[1]

        elif sp_type == 'Species':
            name = 'Species ' + label

        # store in the dictionary
        block = '\n'.join(lines[init_i[i]+1:final_i[i]])
        block = name + '\n' + block
        species_blocks[label].append(block)

    return species_blocks


def dct_species_fragments(species_blocks):
    """ Derive a dictionary of species names and corresponding
        unimol or bimol names.

        :param species_blocks: dictionary with the species blocks
            {name:[frag1 block, frag2 block], name:[unimol block],}
        :type: dict{label: list}
        :return dct_sp_fr: dictionary with species and fragment names
                        {species:[wellname], bimolspecies:[frag1, frag2]}
        :rtype: dct{str: lst}
    """

    dct_sp_fr = dict.fromkeys(species_blocks.keys())

    for label in dct_sp_fr.keys():
        if len(species_blocks[label]) == 1:
            dct_sp_fr[label] = (label,)
        else:
            fragments = ()
            for spc in species_blocks[label]:
                frag = spc.split('\n')[0].split()[1]
                fragments += (frag,)
            dct_sp_fr[label] = fragments

    return dct_sp_fr
