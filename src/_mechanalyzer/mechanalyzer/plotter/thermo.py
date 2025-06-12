""" Plot thermochemical values produced from multiple mechanisms
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy
import math


LINES = ['-', '--', '-.', ':']  # for plot formatting
SORT_DCT = {'h': 1, 'cp': 2, 's': 3, 'g': 4, 'lnq': 5}


def build_plots(algn_spc_therm_dct, spc_dct=None, mech_names=None, 
                sort_method='h', sort_temp=None):
    """ Builds plots of an algn_spc_therm_dct, with one species per page.
        Also plots differences relative to other mechs.

        :param algn_spc_therm_dct: aligned dct with thermo for each mech
        :type algn_spc_therm_dct: dct {spc1: [diff_array_mech1,
            diff_array_mech2, ...], spc2: ...}
        :param spc_dct: spc_dct describing all species in algn_spc_therm_dct
        :type spc_dct: dct {spc1: {spc_info}, spc2: ...}
        :param mech_names: list of mech_names for plot labeling; default is
            'mech1, mech2, ...'
        :type mech_names: list [mech_name1, mech_name2, ...]
        :param sort_method: instructions for sorting; 'h', 'cp', 's', 'g', 'lnq', or None
        :type sort_method: None or str
        :param sort_temp:
        :type sort_temp:
        :return figs: list of MatPlotLib figure objects
        :rtype: list [fig1, fig2, ...]
    """
    # Get the number of mechanisms
    vals = algn_spc_therm_dct.values()
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
            f"""The number of mechs is {num_mechs}, while the length of the
                variable mech_names is {len(mech_names)}.""")

    # Get the algn_spc_diff_dct
    algn_spc_diff_dct = get_algn_spc_diff_dct(algn_spc_therm_dct)
    for spc in algn_spc_diff_dct:
        if algn_spc_diff_dct[spc][1] is None:
            continue
        temps, h_t, cp_t, s_t, g_t, lnq_t = algn_spc_diff_dct[spc][1]
        #sarahs prints
        #for (temp, g) in zip(temps, g_t):
        #    if temp > 650 and temp < 750:
        #        print(spc, temp, g)
    # If indicated, sort the thermo and diff dcts by the differences
    if sort_method:
        algn_spc_diff_dct, algn_spc_therm_dct = sort_by_max_diff(
            algn_spc_diff_dct, algn_spc_therm_dct, sort_method=sort_method,
            sort_temp=sort_temp)

    # Loop over each spc and plot
    figs = []
    for spc, therm_arrays in algn_spc_therm_dct.items():
        diff_arrays = algn_spc_diff_dct[spc]
        if spc_dct is not None:
            if spc_dct.get(spc) is not None:
                smiles = spc_dct.get(spc).get('smiles')
                inchi = spc_dct.get(spc).get('inchi')
                fig, axs = initialize_fig_and_axes(spc, smiles, inchi)
            else:
                fig, axs = initialize_fig_and_axes(spc)
        else:  # if no spc_dct provided, the smiles and inchis will be blank
            fig, axs = initialize_fig_and_axes(spc)
        fig = plot_single_spc(therm_arrays, diff_arrays, fig, axs, mech_names)
        figs.append(fig)

    return figs, algn_spc_therm_dct


def plot_single_spc(therm_arrays, diff_arrays, fig, axs, mech_names):
    """ Plot thermo values of a single species from all mechanisms, as well as
        the differences relative to a reference mechanism

        :param therm_arrays: thermo arrays (T, h, cp, s, g, lnq) for each mech
        :type therm_arrays: list [therm_array_mech1, therm_array_mech2, ...]
        :param diff_arrays: difference arrays (T, h, cp, s, g, lnq) for each mech
        :type diff_arrays: list [diff_array_mech1, diff_array_mech2, ...]
        :param fig: pre-allocated figure object
        :type fig: MatPlotLib figure object
        :param axs: pre-allocated set of 8 axes
        :type axs: list [ax1, ax2, ..., ax8]
        :param mech_names: list of mechanisms names
        :type mech_names: list [mech_name1, mech_name2]
        :return fig: updated figure object with plots added
        :rtype: MatPlotLib figure object
    """
    if len(mech_names) <= 6:
        colors = ['r', 'b', 'g', 'm', 'c', 'y']
    else:  # account for unlikely case where there are more than 6 mechs
        colors = cm.get_cmap('rainbow')(numpy.linspace(0, 1, len(mech_names)))

    diffs_plotted = False
    for mech_idx, therm_array in enumerate(therm_arrays):
        _label = mech_names[mech_idx]
        _color = colors[mech_idx]
        _line = LINES[mech_idx % 4]
        if therm_array is not None:
            temps = therm_array[0]
            # Plot the four thermo values on the even indices plots
            for idx in range(4):
                if idx in (1, 2):  # if s or cp, no divide by 1000
                    therm_ydata = therm_array[idx + 1]
                else:  # otherwise, divide by 1000
                    therm_ydata = therm_array[idx + 1] / 1000
                axs[(idx * 2)].plot(temps, therm_ydata, label=_label,
                                    color=_color, linestyle=_line)

            # Plot the differences if they exist
            if diff_arrays[mech_idx] is not None:
                diff_array = diff_arrays[mech_idx]
                diffs_plotted = True
                # Plot the four diff values on the odd indices plots
                for idx in range(4):
                    if idx in (1, 2):  # if s or cp, no divide by 1000
                        diff_ydata = diff_array[idx + 1]
                    else:  # otherwise, divide by 1000
                        diff_ydata = diff_array[idx + 1] / 1000
                    axs[(idx * 2) + 1].plot(temps, diff_ydata, label=_label, 
                                            color=_color, linestyle=_line)
                    # Set the limits to prevent weird scaling
                    if max(diff_ydata) - min(diff_ydata) < 1:  # if range < 1
                        avg_val = (max(diff_ydata) + min(diff_ydata)) / 2
                        axs[(idx * 2) + 1].set_ylim((
                            avg_val - 1, avg_val + 1))

    # Do some formatting
    for idx in range(4):
        axs[(idx * 2)].legend(fontsize=12, loc='upper right')
    if diffs_plotted:  # prevents annoying warning about legends being empty
        for idx in range(4):
            axs[(idx * 2) + 1].legend(fontsize=12, loc='upper right')
    # If no differences were plotted, plot a white line at 0 to get axis labels
    else:
        for therm_array in therm_arrays:
            if therm_array is not None:  # find an existent therm array
                for idx in range(4):
                    blanks = numpy.zeros_like(therm_array[0])
                    axs[(idx * 2) + 1].plot(therm_array[0], blanks, color='w')
                break  # only need to do once

    return fig


def sort_by_max_diff(algn_spc_diff_dct, algn_spc_therm_dct, sort_method='h',
                     sort_temp=None):
    """ Sort the algn_spc_diff_dct and algn_spc_therm_dct by maximum
        difference between the thermo quantities

        :param algn_spc_diff_dct: aligned dct containing differences in
            thermo quantities relative to a ref mech
        :type algn_spc_diff_dct: dct {spc1: [diff_array_mech1,
            diff_array_mech2, ...], spc2: ...}
        :param algn_spc_therm_dct: aligned dct with thermo for each mech
        :type algn_spc_therm_dct: dct {spc1: [therm_array_mech1,
            therm_array_mech2, ...], spc2: ...}
        :param sort_method: criteria by which to sort; 'h', 'cp', 's', 'g' or 'lnq'
        :type sort_method: str
        :param sort_temp:
        :type sort_temp:
        :return sorted_spc_diff_dct: algn_spc_diff_dct sorted by max diff
        :rtype: dct {spc1: [diff_arr_mech1, diff_arr_mech2, ...], spc2: ...}
        :return sorted_spc_therm_dct: algn_spc_therm_dct sorted by max diff
        :rtype: dct {spc1: [therm_array_mech1, therm_array_mech2, ...],
            spc2: ...}
    """
    # Get the idx of the sort criteria
    sort_idx = SORT_DCT.get(sort_method)
    assert sort_idx is not None, (
        f"sort_method should be 'h', 'cp', 's', 'g' or 'lnq', but is {sort_method}")

    # If a sort_temp was input, get the idx of the sort_temp
    sort_temp_idx = None  # default is None
    if sort_temp is not None:
        # Get the temps array
        temps = None
        for therm_arrays in algn_spc_therm_dct.values():
            for therm_array in therm_arrays:
                if therm_array is not None:
                    temps = therm_array[0]
                    break
        sort_temp_idx = numpy.argmin(abs(temps - sort_temp))

    # Get the maximum differences
    max_diffs = {}
    for spc, diff_arrays in algn_spc_diff_dct.items():
        max_diff = -0.5  # keeps spcs without diff_arrays behind spcs with them
        for diff_array in diff_arrays:
            if diff_array is not None:
                # If a sort_temp was provided, use the corresponding temp_idx
                if sort_temp_idx is not None:
                    if abs(diff_array[sort_idx][sort_temp_idx]) > max_diff:
                        max_diff = abs(diff_array[sort_idx][sort_temp_idx])
                # Otherwise, just look for the maximum among all values
                else:
                    if max(abs(diff_array[sort_idx])) > max_diff:
                        max_diff = max(abs(diff_array[sort_idx]))

        # If no diff_arrays found, check if the spc is missing in some mechs
        if max_diff == -0.5:
            if None not in algn_spc_therm_dct[spc]:
                pass  # if there are no Nones, leave max_ratio
            else:
                for mechidx, therm_array in enumerate(algn_spc_therm_dct[spc]):
                    if therm_array is not None:
                        # Sort by mech_idx
                        max_diff = (mechidx + 1) * -1
                        break
        max_diffs[spc] = max_diff

    # Sort spc keys by max values
    sorted_dct = {}
    for spc, max_diff in sorted(max_diffs.items(),
                                key=lambda item: item[1],
                                reverse=True):
        sorted_dct[spc] = max_diff

    # Reorder the algn_spc_diff_dct and algn_spc_therm_dct
    sorted_spc_diff_dct = {}
    sorted_spc_therm_dct = {}
    for spc in sorted_dct:
        sorted_spc_diff_dct[spc] = algn_spc_diff_dct[spc]
        sorted_spc_therm_dct[spc] = algn_spc_therm_dct[spc]

    return sorted_spc_diff_dct, sorted_spc_therm_dct


def get_algn_spc_diff_dct(algn_spc_therm_dct):
    """ Take an algn_spc_therm_dct and calculate the difference between each
        thermo quantity and the same quantity in a reference mechanism. The
        reference mechanism is the first mechanism in the spc_therm_dct
        that has rate values for that species.

        :param algn_spc_therm_dct: aligned dct with thermo for each mech
        :type algn_spc_therm_dct: dct {spc1: [therm_array_mech1,
            therm_array_mech2, ...], spc2: ...}
        :return: algn_spc_diff_dct: aligned dct containing differences in
            thermo quantities relative to ref mech
        :rtype: dct {spc1: [diff_arr_mech1, diff_arr_mech2, ...], spc2: ...}
    """

    algn_spc_diff_dct = {}
    for spc, therm_arrays in algn_spc_therm_dct.items():
        # Get the ref_therm_array, which is the first non-None array
        # (will fail if any spc has all Nones, but this shouldn't occur)
        for mech_idx, therm_array in enumerate(therm_arrays):
            if therm_array is not None:
                ref_therm_array = therm_array
                ref_mech_idx = mech_idx
                break

        # Get the diff_arrays
        diff_arrays = []
        for mech_idx, therm_array in enumerate(therm_arrays):
            # If on the ref_therm_array or the current therm_array is None, set
            # the diff_array to None
            if mech_idx == ref_mech_idx or therm_array is None:
                diff_array = None
            # Otherwise, calculate the diff_dct
            else:
                diff_array = []
                for idx, quantity in enumerate(therm_array):
                    if idx == 0:  # if on the temps, just store them
                        temps = quantity
                    elif idx == 1:
                        h_t = quantity - ref_therm_array[idx]
                    elif idx == 2:
                        cp_t = quantity - ref_therm_array[idx]
                    elif idx == 3:
                        s_t = quantity - ref_therm_array[idx]
                    elif idx == 4:
                        g_t = quantity - ref_therm_array[idx]
                    elif idx == 5:
                        lnq_t = quantity - ref_therm_array[idx]
                diff_array = (temps, h_t, cp_t, s_t, g_t, lnq_t)
            diff_arrays.append(diff_array)
        algn_spc_diff_dct[spc] = diff_arrays

    return algn_spc_diff_dct


def initialize_fig_and_axes(spc, smiles=None, inchi=None):
    """ Creates the figure and 8 axes for plotting thermo of a species

        :param spc: species name
        :type spc: str
        :param smiles: SMILES string for the species
        :type smiles: str
        :param inchi: InChI string for the species
        :type inchi: str
    """
    fig = plt.figure(figsize=(8.5, 11))
    fig.suptitle(spc, x=0.5, y=0.94, fontsize=20)
    rows = 7
    columns = 2
    grid = plt.GridSpec(rows, columns, wspace=0.3, hspace=0.12)

    # For formatting
    plot_placements = [(0, 2, 0), (2, 3, 0), (4, 6, 0), (6, 7, 0), (0, 2, 1),
                       (2, 3, 1), (4, 6, 1), (6, 7, 1)]
    plot_titles = ('Enthalpy (kcal/mol)', 'Heat capacity (cal/mol-K)',
                   'Entropy (cal/mol-K)', 'Gibbs free energy (kcal/mol)')

    # Create each axis, with some formatting
    axs = []
    for ax_num in range(8):
        # Set the placement and size of the plots
        row_start = plot_placements[ax_num][0]
        row_end = plot_placements[ax_num][1]
        column = plot_placements[ax_num][2]
        axs.append(plt.subplot(grid[row_start:row_end, column]))

        # Apply formatting for the main plots, which are the even indices
        if ax_num % 2 == 0:
            axs[ax_num].set_xticks([])
            axs[ax_num].set_title(plot_titles[int(ax_num/2)])

        # Apply formatting for the residual plots, which are the odd indices
        else:
            axs[ax_num].set_xlabel('Temperature (K)')

    # Add some annotations
    if smiles is None:
        smiles = 'not specified'
    if inchi is None:
        inchi = 'not specified'
    footnotes = f'SMILES: {smiles}\nInChi: {inchi}'
    header = 'Large plots: values\nSmall plots: residuals'
    plt.figtext(0.5, 0.05, footnotes, fontsize=10, va="top", ha="center")
    plt.figtext(0.02, 0.98, header, fontsize=8, va="top", ha="left")

    return fig, axs
