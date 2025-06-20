""" Read the merged wells from the auxiliary file
"""

import numpy
from phydat import phycon
import autoparse.pattern as app
import autoparse.find as apf


def merged_wells(mess_aux_str, pressure, temp):
    """ Parse the auxiliary MESS output file string for all of the groups
        of wells which merged at a given pressure and temperature
        from a Master Equation simulation

        :param mess_aux_str: string for the auxiliary file
        :type mess_aux_str: str
        :rtype: tuple(tuple(str))
    """

    # Convert input pressure to bar for reading
    pressure *= phycon.ATM2BAR

    # Get the lines
    mess_lines = mess_aux_str.splitlines()

    # Parses number of wells and lines they appear on for second loop
    # Done for each set of wells at each T and P
    merged_well_lines = []
    for i, line in enumerate(mess_lines):
        if 'number of species =' in line:
            cond_line = mess_lines[i-1].strip().split()
            well_pressure = float(cond_line[2])
            well_temp = float(cond_line[-2])
            if (
                numpy.isclose(pressure, well_pressure) and
                numpy.isclose(temp, well_temp)
            ):
                nspc = int(line.strip().split()[-1])
                merged_well_lines.append((nspc, i))

    # Find all the lines in the block to grab the wells
    merged_well_lst = ()
    for nspc, line_idx in merged_well_lines:
        for idx in range(nspc):
            well_line = mess_lines[line_idx+1+idx*2]
            well_names = well_line.strip().split()
            if len(well_names) > 1:
                merged_well_lst += (tuple(well_names),)
  
    if not merged_well_lst:
        merged_well_lst = None

    return merged_well_lst


def well_thermal_energy(log_str, well, temp):
    """ Obtain the thermal energy of each well from the output
        of MESS rate calculations.

        Returns the energies in hartrees.

        :param log_str: string of the MESS .log file
        :type log_str: str
        :param well: name of the well to get energy for
        :type well: str
        :param temp: temperature to get the energy for
        :type temp: float
        :rtype: dict[str: float]
    """

    # Loop through file to find the energy block with requested temp
    block_ptt = (
        'MasterEquation::set:  starts' +
        app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
        'MasterEquation::set:  done'
    )
    blocks = apf.all_captures(block_ptt, log_str)

    if blocks is not None:
        ptt = ('Temperature' + app.SPACES + '=' + app.SPACES +
               app.capturing(app.NUMBER) + app.SPACES + 'K')
        for i, block in enumerate(blocks):
            blocktemp = float(apf.first_capture(ptt, block))
            if numpy.isclose(blocktemp, temp, atol=0.01):
                block_idx = i
                break
    else:
        block_idx = None

    # If block with requested temp found, get energies
    ene = None
    if block_idx is not None:
        ptt = (
            app.escape(well) + app.SPACES +
            app.escape('Well:') + app.SPACES +
            'thermal energy =' + app.SPACES +
            app.capturing(app.NUMBER) + app.SPACES +
            app.escape('kcal/mol')
        )
        cap = apf.first_capture(ptt, blocks[block_idx])
        if cap is not None:
            ene = float(cap) * phycon.KCAL2EH
        else:
            # check if user is using old version off MESS, which
            # uses 'average' instead of 'thermal'
            ptt = (
                app.escape(well) + app.SPACES +
                app.escape('Well:') + app.SPACES +
                'average energy =' + app.SPACES +
                app.capturing(app.NUMBER) + app.SPACES +
                app.escape('kcal/mol')
            )
            cap = apf.first_capture(ptt, blocks[block_idx])
            if cap is not None:
                ene = float(cap) * phycon.KCAL2EH
                print(
                    "WARNING: old version of MESS detected from" +
                    "use of 'average' instead of 'thermal' energy")
    return ene
