""" Analyze the input and output file strings of a base MESS calculation
    to determine the parameters required to define a lumped well-extension
    scheme. Use these parameters to write the input file.
"""

import itertools
import numpy
import ioformat
from phydat import phycon
import mess_io.reader

# Main calling function
def well_lumped_input_file(inp_str, out_str, aux_str, log_str,
                           lump_pressure, lump_temp, lump = True):
    """ Run MESS to get the wells and then parse the aux file for wells...
    """

    # Get the requisite information from analyzing the output
    well_enes_dct = well_energies(out_str, log_str, lump_pressure)
    
    well_lump_str = None
    if lump == True and len(well_enes_dct.keys()) > 1:
        # need to have at least 2 wells for merging
        well_lump_str = well_lumping_scheme(
            aux_str, out_str, lump_pressure, lump_temp)
        
    # Write a new string containing the parsed information
    well_extend_str = _format_well_extension_inp(
        inp_str, well_enes_dct, well_lump_str)

    return well_extend_str


# Build the lumping scheme used for the well extension
def well_lumping_scheme(mess_aux_str, out_str, pressure, temp):
    """ Parse lumped wells from aux output; write into string for new input
    """
    # CURRENTLY NOT WORKING 
    # ^I, Sarah, didn't write the comment above
    # but I have implemented a fix so that the the names are correct
    # was this what the comment above was refering to or is 
    # it broken in some other way that I haven't found yet?
    well_lump_lst = mess_io.reader.merged_wells(mess_aux_str, pressure, temp)
    if well_lump_lst is not None:
        lbl_dct = mess_io.reader._label.name_label_dct(out_str)
        well_lump_lst = [
            tuple(lbl_dct.get(lbl, lbl) for lbl in lump_set) for lump_set in well_lump_lst
        ]
        well_lump_str = mess_io.writer.well_lump_scheme(well_lump_lst, separator='&')
    else:
        print(f'No wells are merged at P = {pressure} atm and T = {temp} K')
        well_lump_str = None

    return well_lump_str


# Get the energies for definining the well extension cap
def well_energies(mess_out_str, mess_log_str, pressure):
    """ Obtain the energies for each well at the given pressure.
    """

    mess_temps, _ = mess_io.reader.rates.temperatures(mess_out_str)
    max_run_temp = max(mess_temps)

    # Get the temps where each well exists
    well_enes = {}
    well_rxns = _get_well_reactions(mess_out_str)
    for well, rxn_lst in well_rxns.items():
        print('\n***********************************************\n')
        print(f'Obtaining information for well {well} at P={pressure}')
        max_temp = -1.0
        for rxn in rxn_lst:

            well, prd = rxn[0][0], rxn[1][0]

            # Read the rate constants out of the mess outputs
            ktp_dct = mess_io.reader.rates.ktp_dct(mess_out_str, well, prd)
            rxn_temp = _max_temp_well_exists(ktp_dct, pressure, mess_temps)

            if rxn_temp > max_temp:
                max_temp = rxn_temp
                print(f'- New max temperature for well: {rxn_temp} K')
                print(f'- from reaction {well}->{prd}')

        # Determine if k(T) exist at highest T -> no Well cap exists
        if numpy.isclose(max_temp, max_run_temp):
            max_temp = None
            print('\nMax temperature at highest value from run.',
                  'No cap needed.')
        else:
            print(f'\nMax temperature for energies is {max_temp} K')

        # Read the thermal energy at the max temperature
        # Only put enes if they are positive to be written later
        if max_temp is not None:
            ene = mess_io.reader.well_thermal_energy(mess_log_str, well, max_temp)

            if ene > 0.0:
                well_enes[well] = ene
            else:
                well_enes[well] = None
        else:
            well_enes[well] = None

    # relabel if needed
    lbl_dct = mess_io.reader._label.name_label_dct(mess_log_str)

    if lbl_dct is not None:    
        well_enes_new = {lbl_dct[key]: val for key, val in well_enes.items()}
        well_enes_new = {}
        for key, val in well_enes.items():
            new_key = lbl_dct[key]
            well_enes_new[new_key] = val
    else:
        well_enes_new = well_enes
    
    #filter fake wells to avoid setting wellcaps (reduces their extension)
    for name in well_enes_new.keys():
        if 'FakeW' in name:
            well_enes_new[name] = None
            
    return well_enes_new


def _get_well_reactions(mess_out_str):
    """ For each Well in the output file, Generate a list of all
        reactions where that Well is the reactant
    """

    # Get the well labels from the reactions
    # We assume wells are MESS labels missing a '+' or 'W'
    rxns = mess_io.reader.rates.reactions(mess_out_str)
    rxns = mess_io.reader.rates.filter_reactions(
        rxns, filter_reverse=False)
    
    lbl_dct = mess_io.reader._label.name_label_dct(mess_out_str)
        
    wells = ()
    for rxn in rxns:
        rcts, prds = rxn[0], rxn[1]
        if lbl_dct is None:
            wells += tuple(rct for rct in rcts if '+' not in rct)
            wells += tuple(prd for prd in prds if '+' not in prd)
        else: #renamed 
            wells += tuple(rct for rct in rcts if rct[0] == 'W')
            wells += tuple(prd for prd in prds if prd[0] == 'W')            

    # Remove duplicates from above loop
    wells = tuple(n for i, n in enumerate(wells) if n not in wells[:i])

    # Grab reactions that contains the well as the reactant
    well_rxns = {}
    for well in wells:
        rxn_lst = ()
        for rxn in rxns:
            if (well,) == rxn[0]:
                rxn_lst += (rxn,)
        well_rxns[well] = rxn_lst

    return well_rxns


def _max_temp_well_exists(ktp_dct, pressure, mess_temps):
    """ For a given reaction and pressure, find a max temperature

        Assumes a filtered ktp dct.
    """

    max_temp = None
    for _pressure, kt_lst in ktp_dct.items():
        if _pressure != 'high':
            if numpy.isclose(_pressure, pressure):
                # max_temp = max((t for t in kt_lst[0] if t is not None))
                max_temp = max(kt_lst[0])
                break

    # Default to the lowest temperature run if no rate constants found
    if max_temp is None:
        max_temp = min(mess_temps)
        print(f'\nNo k(T) values found for P = {pressure} atm.')
        print(f'T={max_temp}: minimum of all temps in output')

    return max_temp


# Handlies writing the new string
def _format_well_extension_inp(inp_str, well_enes_dct, well_lump_str):
    """ handles building new input will well lumping/extension info
    """

    # Reinitialize string and uncomment wellext if needed
    if '!WellExtension' in inp_str:
        new_inp_str = inp_str.replace('!WellExtension', 'WellExtension')
    else:
        # add anyways
        new_inp_str = ioformat.add_line(
            string=inp_str, addline='WellExtension',
            searchline='Model', position='before')
        
    # Write string for each of the well enes
    for well, ene in well_enes_dct.items():
        if ene is not None:
            # Find line for where well start, for-loop handle weird format
            for line in inp_str.splitlines():
                if 'Well' in line and well in line:
                    _search = line
                    break       
            _add = f'  WellExtensionCap[kcal/mol]    {ene*phycon.EH2KCAL:.2f}'
            new_inp_str = ioformat.add_line(
                string=new_inp_str, addline=_add,
                searchline=_search, position='after')

    # Write new strings with the lumped input
    well_extend_line = 'ExtensionCorrection    0.6'
    new_inp_str = ioformat.add_line(
        string=new_inp_str, addline=well_extend_line,
        searchline='Model', position='before')

    if well_lump_str is not None:
        well_lump_line = ioformat.indent(well_lump_str, 2)
        new_inp_str = ioformat.add_line(
            string=new_inp_str, addline=well_lump_line,
            searchline='Model', position='after')

    return new_inp_str
