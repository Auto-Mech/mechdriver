""" Processes info in mess format and extracts species attributes
"""

import numpy as np
import pandas as pd
from autoparse import find
from automol import geom
from automol import chi
from automol import form

# MW
MW_dct_elements = {
    'C': 12e-3,
    'N': 14e-3,
    'O': 16e-3,
    'H': 1e-3,
    'S': 32e-3,
    'P': 31e-3,
    'F': 19e-3,
    'Cl': 35.45e-3
}  # kg/mol

def get_info(block, get_ts = True):
    """ Gets the N of degrees of freedom and MW of each species
        :param block: bimol species of which you want the dofs
        :type block: list(str1, str2)
        :param ask_for_ts: build the dof info also for the ts
        :type ask_for_ts: bool
        :return dof_info: dataframe with vibrat/rot degrees of freedom
            and molecular weight
        :rtype: dataframe(index=species, columns=['n_atoms', 'vib dof', 'rot dof', 'mw',
                                                'geometry', 'symmetry', 'freqs', 'hr'])
        'hr' still to implement
    """
    dof_info = pd.DataFrame(index = np.arange(len(block) + 1*get_ts), columns=[
                            'name', 'n_atoms', 'vib dof', 'rot dof', 'mw',
                            'geometry', 'symmetry', 'freqs', 'hr'],
                            dtype = object)
    
    atoms_ts = 0
    hr_dct = {}
    # extract N of dofs and MW
    for i, block_i in enumerate(block):
        info = block_i.splitlines()
        where_name = find.where_in('Species', info)[0]
        where_hind = find.where_in('Hindered', info)
        
        name = info[where_name].strip().split()[1]
        dof_info.loc[i,'name'] = name
        
        if 'Frequencies' in block_i: # more than 1 atom
            where_geom = find.where_in('Geometry', info)[0]
            where_freq = find.where_in('Frequencies', info)[0]
            where_zeroen = find.where_in('ZeroEnergy', info)[0]
            where_rotors = find.where_in('Electronic', info)[0]
            where_symmtot = find.where_in('SymmetryFactor', info)[0]
            where_elecen = find.where_in('Electronic', info)[0]

            num_atoms = int(info[where_geom].strip().split()[1])
            vib_dof = (
                int(info[where_freq].strip().split()[1]) + len(where_hind)
            )
            lines_of_relevantinfo = np.array([where_geom, where_zeroen, where_symmtot, where_elecen])
            lines_of_relevantinfo = np.concatenate((lines_of_relevantinfo, where_hind))
            # save freqs
            # finish reading freqs
            end_freqlines = min(lines_of_relevantinfo[lines_of_relevantinfo > where_freq])
            freqlines = info[where_freq+1:end_freqlines]
            freqs = [frline.strip().split() for frline in freqlines]
            freqsarr = []
            list(map(freqsarr.extend, freqs))
            freqsarr = np.array(freqsarr, dtype=np.float32)
            dof_info.loc[i,'freqs'] = freqsarr
            
            # save hrs
            for hrn, linehr in enumerate(where_hind):
                hr_dct[hrn] = {}
                # start, end
                endhrn = find.where_in('End', info[linehr:])[0]
                hrlines = info[linehr+1:endhrn]
                for nl, line in enumerate(hrlines):
                    key = line.strip().split()[0]
                    if 'Potential' in key:
                        pot = np.array([potline.strip().split() for potline in hrlines[nl+1:]], dtype=np.float32)
                        val = pot.flatten()
                        hr_dct[hrn][key] = val
                        break
                    else:
                        val = [int(value) for value in line.strip().split()[1:]]
                        hr_dct[hrn][key] = val

            if 3*num_atoms - vib_dof == 6:
                rot_dof = 3
            else:
                rot_dof = 2
        # save symm
            dof_info.loc[i,'symmetry'] = float(info[where_symmtot].strip().split()[1])
            
        else: # 1 atom only
            # if 1 atom only: no 'Frequencies', set to 0
            vib_dof = 0
            rot_dof = 0
            num_atoms = 1
            if 'Name' in block_i:
                where_geom = find.where_in('Name', info)[0]
            else:
                where_geom = find.where_in('Mass[amu]', info)[0]
            atoms_array = np.array([info[where_geom].strip().split()[1]])
            
        atoms_ts += num_atoms
        
        # this allows to get 3N-5 or 3N-6 without analyzing the geometry
        dof_info.loc[i,['n_atoms', 'vib dof', 'rot dof', 'symmetry']] = [num_atoms, vib_dof, rot_dof, 1]

        # MW from type of atoms:
        if num_atoms > 1:
            geom_in = where_geom+1
            geom_fin = geom_in+num_atoms
            # freq qui
            geom_array = [(geomline.strip().split()[0], 
                                tuple(np.array(geomline.strip().split()[1:], dtype=np.float32)/0.529)) 
                                for geomline in info[geom_in:geom_fin]]	
            atoms_array = np.array([geomi[0] for geomi in geom_array])    
            dof_info.loc[i,'geometry'] = tuple(geom_array)
            dof_info.loc[i,'mw'] = np.sum(np.array([MW_dct_elements[at]
                                                for at in atoms_array], dtype=float))
            
    # ts info: assume first 2 blocks are 2 reactants of bimol reaction
    # and derive the DOFs of the TS
    if get_ts:
        # bimol rxns
        if len(block) == 2:
            mwts = dof_info.loc[i, 'mw'] + dof_info.loc[i-1, 'mw']
        elif len(block) == 1:
            mwts = dof_info.loc[i, 'mw'] 
        dof_info.loc[i+1,['name', 'n_atoms',
                           'vib dof', 'rot dof', 'mw']] = ['TS', atoms_ts, 
                                                           3*atoms_ts - 7, 3, mwts]
    # reindex
    dof_info = dof_info.set_index('name')

    return dof_info

def get_dof_info_fromspcdct(sp, spc_dct):
    """ Gets the N of degrees of freedom and MW of each species
        :sp: species name
        :type sp: str
        :param spc_dct: species dictionary
        :return dof_info: dictionary with vibrat/rot degrees of freedom
            and molecular weight
        :rtype: dct(['n_atoms', 'vib dof', 'rot dof', 'mw'])
    """
    dof_info = {}
    fml = spc_dct[sp]['fml']
    Nat = form.atom_count(fml)
    dof_info['n_atoms'] = Nat
    dof_info['mw'] = sum(np.array([form.element_count(fml, at) *
                       MW_dct_elements[at] for at in MW_dct_elements.keys()]))
    if Nat == 1:
        dof_info['vib dof'] = 0
        dof_info['rot dof'] = 0
    elif Nat == 2:
        dof_info['vib dof'] = 1
        dof_info['rot dof'] = 2
    else:
        # derive geometry and check if linear
        geom_sp = chi.geometry(spc_dct[sp]['inchi'])
        try:
            ilin = int(geom.is_linear(geom_sp))
        except AssertionError:
            # failed for some reason .. set to 0. check HCO, fails there 
            ilin = 0
        dof_info['rot dof'] = 3 - 1*ilin
        dof_info['vib dof'] = 3*Nat - 6 + 1*ilin
        
    return dof_info




# helper functions for energy partition used in the sorter


