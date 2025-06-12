""" Plot the contents of an algn_rxn_ktp_dct, with comparison
    to other mechanisms as appropriate
"""

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
import numpy
from chemkin_io.writer import _util as writer


LINES = ['-', '--', '-.', ':']  # for plot formatting
K_UNITS_DCT = {1: '(s$^{-1}$)',
               2: '(cm$^3$ mol$^{-1}$ s$^{-1}$)',
               3: '(cm$^6$ mol$^{-2}$ s$^{-1}$)',
               4: '(cm$^9$ mol$^{-3}$ s$^{-1}$)'}


def build_plots(algn_rxn_ktp_dct, mech_names=None, ratio_sort=False):
    """ Build plots of an algn_rxn_ktp_dct, with one reaction per page. Also calculate ratios
        relative to other mechs and plot the ratios. Output a PDF.

        :param algn_rxn_ktp_dct: aligned dct containing rates for each mech
        :type algn_rxn_ktp_dct: dct {rxn1: [ktp_dct_mech1, ktp_dct_mech2, ...], rxn2: ...}
        :param mech_names: list of mech_names for plot labeling; default is 'mech1, mech2, ...'
        :type mech_names: list [mech_name1, mech_name2, ...]
        :param ratio_sort: whether or not to sort plots by the max value of the ratio
        :type ratio_sort: Bool
        :return figs: list of MatPlotLib figure objects
        :rtype: list [fig1, fig2, ...]
    """
    # Get the number of mechanisms
    vals = algn_rxn_ktp_dct.values()
    val_iter = iter(vals)
    first_val = next(val_iter)
    num_mechs = len(first_val)

    # Check/fix the mech_names input
    if mech_names is None:
        mech_names = []
        for i in range(num_mechs):
            mech_names.append('mech' + str(i + 1))
    else:
        assert num_mechs == len(mech_names), (
            f"""The number of mechs is {num_mechs}, while the length of the variable
                   mech_names is {len(mech_names)}."""
        )

    # Get the algn_rxn_ratio_dct
    algn_rxn_ratio_dct = get_algn_rxn_ratio_dct(algn_rxn_ktp_dct)

    # If indicated, sort the ktp and ratio dcts by the ratio value
    if ratio_sort:
        algn_rxn_ratio_dct, algn_rxn_ktp_dct = sort_by_max_ratio(
            algn_rxn_ratio_dct, algn_rxn_ktp_dct
        )

    # Loop over each rxn and plot
    figs = []
    pressures = get_pressures(algn_rxn_ktp_dct)
    format_dct = get_format_dct(pressures)  # defines color and label for each pressure
    for rxn, ktp_dcts in algn_rxn_ktp_dct.items():
        ratio_dcts = algn_rxn_ratio_dct[rxn]
        molecularity = get_molecularity(rxn)
        fig, axs = build_fig_and_axs(molecularity, ratio_dcts, mech_names)
        fig = plot_single_rxn(rxn, ktp_dcts, ratio_dcts, fig, axs, mech_names, format_dct)
        figs.append(fig)

    return figs


def plot_single_rxn(rxn, ktp_dcts, ratio_dcts, fig, axs, mech_names, format_dct):
    """ Plot a single reaction's k(T,P) values from all mechanisms and the ratio values relative to
        a reference mechanism (the reference mechanism is the first mechanism in the ktp_dct that
        has rate values for that pressure).

        :param rxn: rxn tuple describing the reaction
        :type rxn: tuple (rcts, prds, third_bods)
        :param ktp_dcts: list of ktp_dcts, one for each mechanism (some may be None)
        :type ktp_dcts: list [ktp_dct_mech1, ktp_dct_mech2, ...]
        :param ratio_dcts: list of ratio_dcts, one for each mechanism (some may be None)
        :type ratio_dcts: list [ratio_dct_mech1, ratio_dct_mech1, ...]
        :param fig: pre-allocated figure object
        :type fig: MatPlotLib figure object
        :param axs:
        :type axs: list [ax1, ax2]
        :param mech_names:
        :type mech_names: list [mech_name1, mech_name2]
        :param format_dct: dct containing color and label for each pressure
        :type: dct {pressure1: (color1, label1), pressure2: ...}
    """

    ratios_plotted = False
    for mech_idx, ktp_dct in enumerate(ktp_dcts):
        if ktp_dct is not None:
            for pressure, (temps, kts) in ktp_dct.items():
                (_color, _label) = format_dct[pressure]
                _label += ', ' + mech_names[mech_idx]

                # Plot the rate constants
                axs[0].plot(1000 / temps,
                            numpy.log10(kts), label=_label,
                            color=_color,
                            linestyle=LINES[mech_idx])

                # Plot the ratios if they exist
                if ratio_dcts[mech_idx] is not None:
                    # if ratio_dcts[mech_idx][pressure] is not None:
                    if pressure in ratio_dcts[mech_idx].keys():
                        (_, ratios) = ratio_dcts[mech_idx][pressure]
                        ratios_plotted = True
                        axs[1].plot(1000 / temps,
                                    numpy.log10(ratios), label=_label,
                                    color=_color,
                                    linestyle=LINES[mech_idx])
            # Check for the 'max_to_high' case
            # Grab set of temps to calculate ratio
            # BELOW CODE ASSUMES ALL TEMP RANGES IN KTP DCT THE SAME
            ratio_temps = tuple(ktp_dct.values())[0][1]
            if ratio_dcts[mech_idx] is not None:
                if 'max_to_high' in ratio_dcts[mech_idx].keys():
                    (_, ratios) = ratio_dcts[mech_idx]['max_to_high']
                    ratios_plotted = True
                    _color = 'k'
                    _label = 'max to P-indep, ' + mech_names[mech_idx]
                    axs[1].plot(1000/ratio_temps, numpy.log10(ratios),
                                label=_label, color=_color,
                                linestyle=LINES[mech_idx])

    # Do some formatting
    axs[0].legend(fontsize=12, loc='upper right')
    if ratios_plotted:
        axs[1].legend(fontsize=12, loc='upper right')
    rxn_name_formatted = writer.format_rxn_name(rxn)
    fig.suptitle(rxn_name_formatted, x=0.5, y=0.94, fontsize=20)

    return fig


def sort_by_max_ratio(algn_rxn_ratio_dct, algn_rxn_ktp_dct):
    """ Sort the algn_rxn_ratio_dct and algn_rxn_ktp_dct by maximum ratio value

        :param algn_rxn_ratio_dct: aligned dct containing ratios of rates relative to a ref mech
        :type algn_rxn_ratio_dct: dct {rxn1: [ratio_dct_mech1, ratio_dct_mech2, ...], rxn2: ...}
        :param algn_rxn_ktp_dct: aligned dct containing rates for each mech
        :type algn_rxn_ktp_dct: dct {rxn1: [ktp_dct_mech1, ktp_dct_mech2, ...], rxn2: ...}
        :return sorted_rxn_ratio_dct: algn_rxn_ratio_dct sorted by maximum ratio
        :rtype: dct {rxn1: [ratio_dct_mech1, ratio_dct_mech2, ...], rxn2: ...}
        :return sorted_rxn_ktp_dct: algn_rxn_ktp_dct sorted by maximum ratio
        :rtype: dct {rxn1: [ktp_dct_mech1, ktp_dct_mech2, ...], rxn2: ...}
    """
    # Get the maximum ratios
    max_ratios = {}
    for rxn, ratio_dcts in algn_rxn_ratio_dct.items():
        max_ratio = -0.5  # this value keeps rxns without ratio_dcts behind rxns with them
        for ratio_dct in ratio_dcts:
            if ratio_dct is not None:
                for (_, ratios) in ratio_dct.values():
                    if max(abs(numpy.log10(ratios))) > max_ratio:  # abs considers negative ratios
                        max_ratio = max(abs(numpy.log10(ratios)))
        # If no ratio_dcts were found, check if the rxn is missing in some mechs; order accordingly
        if max_ratio == -0.5:
            if None not in algn_rxn_ktp_dct[rxn]:  # if there are no Nones, leave max_ratio
                pass
            else:
                for mech_idx, ktp_dct in enumerate(algn_rxn_ktp_dct[rxn]):
                    if ktp_dct is not None:
                        max_ratio = (mech_idx + 1) * -1  # find the first mech with non-None
                        break
        max_ratios[rxn] = max_ratio

    # Sort rxn keys by max values
    sorted_dct = {}
    for rxn, max_ratio in sorted(max_ratios.items(),
                                 key=lambda item: item[1],
                                 reverse=True):
        sorted_dct[rxn] = max_ratio
    # Reorder the algn_rxn_ratio_dct and algn_rxn_ktp_dct
    sorted_rxn_ratio_dct = {}
    sorted_rxn_ktp_dct = {}
    for rxn in sorted_dct:
        sorted_rxn_ratio_dct[rxn] = algn_rxn_ratio_dct[rxn]
        sorted_rxn_ktp_dct[rxn] = algn_rxn_ktp_dct[rxn]

    return sorted_rxn_ratio_dct, sorted_rxn_ktp_dct


def get_algn_rxn_ratio_dct(algn_rxn_ktp_dct):
    """ Take an algn_rxn_ktp_dct and calculate the ratio of each rate relative to a reference
        mechanism; the output is an algn_rxn_ratio_dct. The reference mechanism is the first
        mechanism in the ktp_dct that has rate values for that pressure.

        :param algn_rxn_ktp_dct: aligned dct containing rates for each mech
        :type algn_rxn_ktp_dct: dct {rxn1: [ktp_dct_mech1, ktp_dct_mech2, ...], rxn2: ...}
        :return: algn_rxn_ratio_dct: aligned dct containing ratios of rates relative to ref mech
        :rtype: dct {rxn1: [ratio_dct_mech1, ratio_dct_mech2, ...], rxn2: ...}
    """
    algn_rxn_ratio_dct = {}
    for rxn, ktp_dcts in algn_rxn_ktp_dct.items():
        ref_ktp_dct = ktp_dcts[0]
        ratio_dcts = []
        for mech_idx, ktp_dct in enumerate(ktp_dcts):
            # If (1) on first ktp_dct, (2) ref_ktp_dct is None, or (3) current
            # ktp_dct is None, set the ratio_dct to None
            if mech_idx == 0 or ref_ktp_dct is None or ktp_dct is None:
                ratio_dct = None
            # Otherwise, calculate the ratio_dct
            else:
                ratio_dct = {}
                for pressure, (temps, kts) in ktp_dct.items():
                    # If pressure defined in ref ktp_dct, calculate the ratio
                    if pressure in ref_ktp_dct.keys():
                        _, ref_kts = ref_ktp_dct[pressure]
                        ratios = kts / ref_kts
                        ratio_dct[pressure] = (temps, ratios)
                # If the ratio_dct is still empty, one has 'high' and one has
                # numerical pressures (e.g., PLOG and Troe); calculate
                # using the highest numerical pressure

                # Grab set of temps to calculate ratio
                # BELOW CODE ASSUMES ALL TEMP RANGES IN KTP DCT THE SAME
                ratio_temps = tuple(ktp_dct.values())[0][1]
                if ratio_dct == {}:
                    if 'high' in ref_ktp_dct.keys():
                        _, ref_kts = ref_ktp_dct['high']
                        max_pressure = max(ktp_dct.keys())
                        _, kts = ktp_dct[max_pressure]
                        ratios = kts / ref_kts
                        ratio_dct['max_to_high'] = (ratio_temps, ratios)
                    elif 'high' in ktp_dct.keys():
                        max_pressure = max(ref_ktp_dct.keys())
                        _, ref_kts = ref_ktp_dct[max_pressure]
                        _, kts = ktp_dct['high']
                        ratios = kts / ref_kts
                        ratio_dct['max_to_high'] = (ratio_temps, ratios)
                # Catch case when ratio_dct is still empty
                if ratio_dct == {}:
                    ratio_dct = None
            ratio_dcts.append(ratio_dct)  # store
        algn_rxn_ratio_dct[rxn] = ratio_dcts  # store

    return algn_rxn_ratio_dct


def get_pressures(algn_rxn_ktp_dct):
    """ Get all the unique  pressure values in an algn_rxn_ktp_dct
        ('high' is excluded).

        :param algn_rxn_ktp_dct: aligned dct containing rates for each mech
        :type algn_rxn_ktp_dct: dct {rxn1: [ratio_dct_mech1, ratio_dct_mech2, ...], rxn2: ...}
        :return pressures: pressures included in the rxn_ktp_dct (not including 'high')
        :rtype: list [pressure1, pressure2, ...]
    """
    pressures = []
    for ktp_dcts in algn_rxn_ktp_dct.values():
        for ktp_dct in ktp_dcts:
            if ktp_dct is not None:
                for pressure in ktp_dct.keys():
                    pressures.append(pressure)
    # Get the unique pressures
    pressures = list(set(pressures))
    # Remove the 'high' entry
    if 'high' in pressures:
        pressures.remove('high')

    return pressures


def get_molecularity(rxn):
    """ Get the molecularity of a reaction

        :param rxn: rxn tuple describing the reaction
        :type rxn: tuple (rcts, prds, third_bods)
        :return molecularity: molecularity of reaction; 1=unimol, 2=bimol, etc.
        :rtype: int
    """
    rcts, _, third_bods = rxn
    num_rcts = len(rcts)
    third_bod = third_bods[0]
    if third_bod is not None:
        # if the third body does not count toward molecularity
        if third_bod[0:2] == '(+':
            molecularity = num_rcts
        else:  # if the third body counts as another species
            molecularity = num_rcts + 1
    # If no third body
    else:
        molecularity = num_rcts

    return molecularity


def get_format_dct(pressures):
    """ Set up the formatting dictionary that describes colors and labels

        :param pressures: pressures at which calculations were performed
        :type pressures: list
        :return format_dct: dct containing color and label for each pressure
        :rtype: dct {pressure1: (color1, label1), pressure2: ...}
    """
    if len(pressures) <= 6:
        colors = ['r', 'b', 'g', 'm', 'c', 'y']
    else:  # account for unlikely case where there are more than 6 pressures
        colors = cm.get_cmap('rainbow')(numpy.linspace(0, 1, len(pressures)))
    format_dct = {}
    for idx, pressure in enumerate(pressures):
        format_dct[pressure] = (colors[idx % 6], str(pressure) + ' atm')
    format_dct['high'] = ('k', 'P-indep')

    return format_dct


def build_fig_and_axs(molecularity, ratio_dcts, mech_names):
    """ Build a figure with two axes

        :param molecularity: molecularity of a reaction
        :type molecularity: int
        :param ratio_dcts: ratio_dcts for this reaction
        :param mech_names: short names for mechs
        :return:
    """
    # Set up the y labels for the rate plot and the ratios plot
    k_units = K_UNITS_DCT[molecularity]
    k_label = 'log$_{10}$$k$ ' + k_units
    ratio_label = ''
    for mech_idx, ratio_dct in enumerate(ratio_dcts):
        if ratio_dct is not None:
            ref_mech_name = mech_names[mech_idx-1]
            ratio_label = ('log$_{10}$ of $k$ ratio relative to' +
                           f' {ref_mech_name}')
            break

    # Set up the figure and axes
    fig = plt.figure(figsize=(8.5, 11))
    axs = []
    axs.append(fig.add_subplot(211))
    plt.xlabel('1000/$T$ (K$^{-1}$)', fontsize=14)
    plt.ylabel(k_label, fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    axs[0].yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
    axs.append(fig.add_subplot(212))
    plt.xlabel('1000/$T$ (K$^{-1}$)', fontsize=14)
    plt.ylabel(ratio_label, fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    axs[1].yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))

    return fig, axs


def plot_k_vs_p(ktp_dct):
    """ Makes a plot of k vs P for a single reaction
    """

    fig = plt.figure(figsize=(8.5, 11))
    axs = []
    axs.append(fig.add_subplot(211))
    plt.xlabel('log$_{10}$(P) (atm)', fontsize=14)
    plt.ylabel('log$_{10}$($k$/$k_{inf}$)', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    axs[0].yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
    # Get the k_of_p numpy array
    k_high = get_k_high(ktp_dct)
    kpt, temps, pressures = get_kpt(ktp_dct)
    # Plot a line for every temperature in the array
    for temp_idx, temp in enumerate(temps):
        plt.plot(numpy.log10(pressures), numpy.log10(kpt[temp_idx] / k_high[temp_idx]), 
                 label=f'{temp:.0f} K')
    axs[0].legend(fontsize=12, loc='upper left')
    
    return [fig]
    

def get_k_high(ktp_dct):
    """ Gets the HPL values as a function of T

        numpy array of shape (num_temps,)
    """ 

    k_high = ktp_dct['high'][1]

    return k_high


def get_kpt(ktp_dct):
    """ Gets an array describing k as a function of P for multiple temps
        Accompanies the plot_k_vs_p function above

        :return kpt: rate constants as function of P for multiple temps
        :type kpt: numpy array of shape (num_temps, num_pressures)
    """

    def get_temps_pressures(ktp_dct):
        """ Reads a ktp_dct and gets the list of pressure and corresponding list of
            temperature arrays
    
            :param ktp_dct: k(T,P) values
            :type ktp_dct: dict {pressure: (temps, kts)}
            :return temps_lst: list of temperature arrays at each pressure (K) 
            :rtype: list [numpy.ndarray1, numpy.ndarray2, ...]
            :return pressures: list of pressures (atm)
            :rtype: list
        """ 
    
        temps_lst = []
        pressures = []
        if 'high' in ktp_dct:
            ktp_dct.pop('high')
        for pressure, (temps, _) in ktp_dct.items():
            temps_lst.append(temps)
            pressures.append(pressure)
    
        return temps_lst, pressures
    
    # Get the rates as functions of pressure
    temps_lst, pressures = get_temps_pressures(ktp_dct)
    temps = temps_lst[0]
    kpt = []
    for temp_idx, temp in enumerate(temps):
        k_of_p = []
        for pressure in pressures:
            pressure
            k_of_p.append(ktp_dct[pressure][1][temp_idx])
        
        kpt.append(k_of_p)

    kpt = numpy.array(kpt)    

    return kpt, temps, pressures

