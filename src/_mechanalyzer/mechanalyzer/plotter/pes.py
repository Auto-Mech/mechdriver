"""
Calls MatPlotLib functionality to create the plot
"""

import numpy
import matplotlib
from matplotlib import pyplot as plt
import igraph
import networkx as nx
# from pyvis.network import Network
import automol.util


# Set plotting options
COLORS = ['k', 'b', 'r', 'g', 'm', 'y', 'c', '#ff9333']
LINESTYLES = ['-', '--', '-.']
MARKERS = ['.', 'o', 's']

# Set various labels for plotting
FIG_TITLE = 'Potential Energy surface'
FONT = {'size': 14}
matplotlib.rc('font', **FONT)

# Parameters for axes set-up
X_AXIS_RIGHT_EXTEND = 0.75
Y_AXIS_TOP_EXTEND = 2.5
Y_AXIS_BOT_EXTEND = 2.5
Y_AXIS_LABEL = 'Relative Enthalpy at 0 K (kcal/mol)'

# Parameters for name and energy labels
NAME_VSHIFT_SCALEF = 0.05
NAME_FONT_SIZE = 18
NAME_LATEX_FORMAT = 'off'
ENE_VSHIFT_SCALEF = 0.05
ENE_FONT_SIZE = 16

# Parameters to draw lines for each species
SPC_LINE_SPACING = 1.15
SPC_LINE_BASE = 0.25
SPC_LINE_LEN = 0.50
SPC_LINE_WIDTH = 4.5
SPC_LINE_COLOR = 'k'
CONN_LINE_WIDTH = 1.5
CONN_LINE_COLOR = 'k'

# Output File Parameters
PLOT_EXT = 'pdf'
PLOT_NAME = 'surface'
PLOT_WIDTH = 16
PLOT_HEIGHT = 9
PLOT_DPI = 1000


# MAKE THE PLOT
def build(ene_dct, conn_lst):
    """ Make a plot of the PES

        :m ene_dct: relative energies for each species
        :type ene_dct: dict[name: energy]
        :m conn_lst: list of all the connections of the PES
        :type conn_lst: lst(str)
    """

    # Generate the coordinates
    spc_coord_dct = _format_coords(ene_dct)

    # Calculate useful information needed for plot
    max_ene, min_ene, spc_cnt = _ranges(ene_dct)
    name_vshift = _calc_vshifts(max_ene, min_ene)
    _, y_axis_tlim, y_axis_blim = _calc_axis_limits(
        max_ene, min_ene, spc_cnt)

    # Build the plot object
    fig, axes = _build_figure()
    _build_axes(axes, y_axis_tlim, y_axis_blim)
    _plt_species_lines(axes, spc_coord_dct, name_vshift)
    _plt_connecting_lines(axes, spc_coord_dct, conn_lst)

    # Create the image
    _create_img(fig)


# FORMAT THE DATA TO PLACE OBJECTS ONTO THE PLOT #
def _format_coords(ene_dct,
                   spc_line_spacing=SPC_LINE_SPACING,
                   spc_line_base=SPC_LINE_BASE,
                   spc_line_len=SPC_LINE_LEN):
    """ Generates the coordinates for each of the species in the list.
        Each species is plotted as a horizontal line with some finite length.

        :m ene_dct: relative energies for each species
        :type ene_dct: dict[name: energy]
        :m spc_line_spacing: distance between each line
        :type spc_line_spacing: float
        :m spc_line_base: minimal distance between each horizontal line
        :type spc_line_base: float
        :m spc_line_len: length of the each horizontal line
        :type spc_line_len: float
        :return spc_coord_dct: coordinates for each species
        :rtype: dict[name: ((x1, x2), y)]
    """

    spc_coord_dct = {}
    for idx, (name, ene) in enumerate(ene_dct.items()):
        xcoord1 = spc_line_spacing * idx + spc_line_base
        xcoord2 = xcoord1 + spc_line_len
        spc_coord_dct[name] = ((xcoord1, xcoord2), ene)

    return spc_coord_dct


# ADD OBJECTS TO THE PES PLOT #
def _build_figure(width=PLOT_WIDTH, height=PLOT_HEIGHT):
    """ Initialize the size and format of the plot figure.

        :m width: width of the output file image (inches)
        :type width: float
        :m height: height of the output file image (inches)
        :type height: float
        :return: figure object for single page of plots
        :rtype: matplotlib.pyplot object
        :return: axes object for single page of plots
        :rtype: matplotlib.pyplot object
    """

    # Initialize plot objects
    fig, axes = plt.subplots(
        nrows=1, ncols=1, figsize=(width, height))

    # Set various plot options
    fig.tight_layout()
    fig.subplots_adjust(left=0.075,
                        top=0.920, bottom=0.075,
                        wspace=0.2, hspace=0.175)

    return fig, axes


def _build_axes(axes_obj, y_axis_tlim, y_axis_blim, tick_intvl=5.0):
    """ Set various things for the axes
    """

    # Set attributes of the y-axis
    axes_obj.set_ylabel(
        'Relative Enthalpy at 0 K (kcal/mol)',
        fontsize='24'
    )

    y_ticks_range = numpy.arange(
        y_axis_blim, y_axis_tlim+tick_intvl, tick_intvl)
    print('tick range', y_ticks_range)
    axes_obj.set_ylim(bottom=y_axis_blim, top=y_axis_tlim)
    axes_obj.set_yticks(y_ticks_range)
    # axes_obj.yaxis.set_major_formatter(
    #     ticker.FormatStrFormatter('%0.2f'))
    axes_obj.set_yticklabels(
        y_ticks_range, fontsize='18'
    )

    # Set attributes of the x-axis
    axes_obj.xaxis.set_ticks([])

    # Remove all of the axis frames except the left one
    axes_obj.spines['bottom'].set_visible(False)
    axes_obj.spines['top'].set_visible(False)
    axes_obj.spines['right'].set_visible(False)


def _plt_species_lines(axes_obj, spc_coord_dct, name_vshift,
                       name_font_size=NAME_FONT_SIZE,
                       species_line_color=SPC_LINE_COLOR,
                       species_line_width=SPC_LINE_WIDTH):
    """ Plot the lines that correspond to the molecular species. Each line
        is a horizontal line.

        :m spc_coord_dct: coordinates of each species line
        :type spc_coord_dct: dict[name: ((x1, x2), y)]
        :m species_line_color: color of the species line
        :type species_line_color: str
        :m species_line_width: (vertical) width of the species line
        :type species_line_width: float
        :m name_vshift: vertical distance that text appears over line
        :type name_vshift: float
        :m name_font_size: size of the font of the name label
        :type name_font_size: float
    """

    for spc_name, ((xp1, xp2), ene) in spc_coord_dct.items():

        # Plot the line for the species
        axes_obj.plot([xp1, xp2], [ene, ene],
                      color=species_line_color,
                      lw=species_line_width, linestyle='-')

        # Plot the label for the line that has the species name
        label_x, label_y = _position_text(xp2, xp1, ene, name_vshift)
        axes_obj.text(label_x, label_y, spc_name,
                      weight='bold', horizontalalignment='center',
                      fontsize=name_font_size, color='black')


def _plt_connecting_lines(axes_obj, spc_coord_dct, conn_lst,
                          color=CONN_LINE_COLOR,
                          line_width=CONN_LINE_WIDTH):
    """ Plot the lines that connect each of the horizontal species lines that
        show how the species are connected.
    """

    for spc1, spc2 in conn_lst:
        (_, xcoord2), ycoord1 = spc_coord_dct[spc1]
        (xcoord1, _), ycoord2 = spc_coord_dct[spc2]

        axes_obj.plot([xcoord2, xcoord1], [ycoord1, ycoord2],
                      color=color, lw=line_width, linestyle='--')


def _create_img(fig,
                width=PLOT_WIDTH, height=PLOT_HEIGHT,
                name=PLOT_NAME, ext=PLOT_EXT, dpi=PLOT_DPI):
    """ Creates the outgoing plot in the format requested by the user.

        :m width: width of the output file image (inches)
        :type width: float
        :m height: height of the output file image (inches)
        :type height: float
        :m name: base name of the output file image
        :type name: str
        :m name: file extension of the output file image
        :type name: str
        :m dpi: dots-per-inch (resolution) of the output file image
        :type dpi: float
    """

    fig_name = f'{name}.{ext}'.format(name, ext)
    fig.set_size_inches(width, height)
    fig.savefig(fig_name, dpi=dpi)
    plt.close(fig)


# HELPER FUNCTIONS
def _calc_vshifts(max_ene, min_ene,
                  name_vshift_scalef=NAME_VSHIFT_SCALEF):
    """ Using input from user, determine meters for plot formatting.
    """

    energy_range = max_ene - min_ene

    name_vshift = name_vshift_scalef * energy_range
    # ene_vshift = ene_vshift_scalef * energy_range

    return name_vshift
    # return name_vshift, ene_vshift


def _calc_axis_limits(max_ene, min_ene, spc_cnt,
                      x_axis_right_extend=X_AXIS_RIGHT_EXTEND,
                      y_axis_top_extend=Y_AXIS_TOP_EXTEND,
                      y_axis_bot_extend=Y_AXIS_BOT_EXTEND):
    """ Using input from user, determine meters for plot formatting. """

    x_axis_rlim = spc_cnt + x_axis_right_extend
    y_axis_tlim = max_ene + y_axis_top_extend
    y_axis_blim = min_ene - y_axis_bot_extend

    return x_axis_rlim, y_axis_tlim, y_axis_blim


def _ranges(ene_dct):
    """ calc limits of the energy
    """
    enes = ene_dct.values()
    return max(enes), min(enes), len(enes)


def _position_text(xp2, xp1, ene, name_vshift):
    """ Set the position of the text
    """

    # Set the positions of the text
    label_x = (xp2 + xp1) / 2.0
    # label_x = ((xp2 + xp1) / 2.0) - ((xp2 - xp1) * 0.1)  # old
    label_y = ene + name_vshift

    return label_x, label_y


# if PlotParameter.name_latex_format == 'on':
#     for k in range(0, Molecule.species_count):
#         tmp = '$' + name[k] + '$'
#         name[k] = tmp


# Sorting functions
def resort_names(ene_dct, conn_lst):
    """ Discern a new sorting for the names of the ene dct
        so that the plotter can make a plot that is readable
    """

    def _connected_partners(spc, conn_lst):
        """ Get names of all PES parts that spc is connected to. """
        conns = []
        for conn in conn_lst:
            if spc in conn:
                lconn, rconn = conn
                partner = lconn if spc != lconn else rconn
                conns.append(partner)
        return conns

    def _remove_listed(name_lst, plot_lst):
        """ Return name_lst entries not already accounted for in plot_lst """
        return [name for name in name_lst if name not in plot_lst]

    def _unchecked_wells(plot_lst, conn_lst, min_well, side):
        """ Find any leftward wells in plot list whose connections have not
            been checked.
        """

        chk_lst = plot_lst if side == 'left' else reversed(plot_lst)

        _wells = []
        for name in chk_lst:
            # Stop at the minimum-well to avoid going to opposite side
            if name == min_well:
                break
            if 'W' in name:
                _wells.append(name)  # wrong need check

        _unchked_wells = []
        for _well in _wells:
            conn_spc = _connected_partners(_well, conn_lst)
            conn_spc = _remove_listed(conn_spc, plot_lst)
            if conn_spc:
                _unchked_wells.append(_well)

        return _unchked_wells

    def _add_spc(well, plot_lst, side):
        """ Add barriers and then products to side
        """
        print(side)
        conn_spc = _connected_partners(well, conn_lst)
        conn_spc = _remove_listed(conn_spc, plot_lst)
        if conn_spc:
            for spc in conn_spc:
                # Add barrier to list
                if side == 'left':
                    plot_lst.insert(0, spc)
                else:
                    plot_lst.append(spc)
                # Add products (wells or products)
                conn_spc2 = _connected_partners(spc, conn_lst)
                conn_spc2 = _remove_listed(conn_spc2, plot_lst)
                for spc2 in conn_spc2:
                    if side == 'left':
                        plot_lst.insert(0, spc2)
                    else:
                        plot_lst.append(spc2)

        return plot_lst

    # Initialize list that will have the PES components in plot order
    plot_names = []

    # Find the deepest well
    min_ene, min_well = 10000.0, None
    for name, ene in ene_dct.items():
        if 'W' in name and ene < min_ene:
            min_ene = ene
            min_well = name
    plot_names.append(min_well)

    # Find the barriers the min-well is connected, add to list
    # conn_spc = _connected_partners(min_well, conn_lst)
    # for i, spc in enumerate(conn_spc):
    #     if (i+1) % 2 == 1:
    #         plot_names.insert(0, spc)
    #     else:
    #         plot_names.append(spc)
    # print('plot names 1', plot_names)

    # Now build out the left and right
    # left_well = None
    # for name in plot_names:
    #     if 'W' in name:
    #         left_well = name
    # right_well = None
    # for name in reversed(plot_names):
    #     if 'W' in name:
    #         right_well = name

    # Build out left side of list
    left_well = min_well
    left_finished = False
    while not left_finished:
        plot_names = _add_spc(left_well, plot_names, 'left')
        unchecked_names = _unchecked_wells(
            plot_names, conn_lst, min_well, 'left')
        if unchecked_names:
            left_well = unchecked_names[0]
        else:
            left_finished = True

    # Build out right side of list
    right_well = min_well
    right_finished = False
    while not right_finished:
        plot_names = _add_spc(right_well, plot_names, 'right')
        unchecked_names = _unchecked_wells(
            plot_names, conn_lst, min_well, 'right')
        if unchecked_names:
            right_well = unchecked_names[0]
        else:
            right_finished = True

    # Build a new ene_dct with the ordered plot names
    ord_ene_dct = {name: ene_dct[name] for name in plot_names}

    print(plot_names)
    print(ord_ene_dct)
    return ord_ene_dct


# Attempt 2 using graphs #
# HAS PYLINT ISSUES WITH THE DRAW FUNCTION
# class GraphArtist(Artist):
#     """Matplotlib artist class that draws igraph graphs.
#
#     Only Cairo-based backends are supported.
#     """
#
#     def __init__(self, graph, bbox, palette=None, *args, **kwds):
#         """Constructs a graph artist that draws the given graph within
#         the given bounding box.
#
#         `graph` must be an instance of `igraph.Graph`.
#         `bbox` must either be an instance of `igraph.drawing.BoundingBox`
#         or a 4-tuple (`left`, `top`, `width`, `height`). The tuple
#         will be passed on to the constructor of `BoundingBox`.
#         `palette` is an igraph palette that is used to transform
#         numeric color IDs to RGB values. If `None`, a default grayscale
#         palette is used from igraph.
#
#         All the remaining positional and keyword arguments are passed
#         on intact to `igraph.Graph.__plot__`.
#         """
#         Artist.__init__(self)
#
#         if not isinstance(graph, Graph):
#             raise TypeError(f'expected igraph.Graph, got {type(graph)}')
#
#         self.graph = graph
#         self.palette = palette or palettes["gray"]
#         self.bbox = BoundingBox(bbox)
#         self.args = args
#         self.kwds = kwds
#
#     def draw(self, renderer):
#         if not isinstance(renderer, RendererCairo):
#             raise TypeError('plotting is supported only on Cairo backends')
#         self.graph.__plot__(renderer.gc.ctx, self.bbox, self.palette,
#                             *self.args, **self.kwds)
#
#
def pes_graph(conn_lst, ene_dct=None, label_dct=None, file_name='surface.pdf'):
    """ Make an igraph argument
    """

    ene_dct, conn_lst, rgts_lst = _format(ene_dct, conn_lst, label_dct)

    # Make an igraph object
    pes_gra = igraph.Graph()

    # Add vertices and set the name and ene attributes
    pes_gra.add_vertices(len(rgts_lst))
    pes_gra.vs["name"] = rgts_lst
    pes_gra.vs["label"] = rgts_lst
    # pes_gra.vs["energy"] = list(ene_dct.values())

    # Write the conn_lst in terms of indices to add the edges
    name_idx_dct = {name: idx for idx, name in enumerate(rgts_lst)}
    idx_conn_lst = ()
    for conn in conn_lst:
        idx_conn_lst += (tuple(name_idx_dct[x] for x in conn),)
    pes_gra.add_edges(idx_conn_lst)

    # Make a plot with some layout
    fig, axes = plt.subplots(
        nrows=1, ncols=1, figsize=(16, 9))
    visual_style = {
        "vertex_size": 30,
        "vertex_label_size": 22,
        # "vertex_color": [color_dict[gender] for gender in g.vs["gender"]],
        "vertex_label": pes_gra.vs["name"],
        # "edge_width": [1+2*int(is_form) for is_form in g.es["is_form"]],
        "layout": "kamada_kawai",
        "bbox": (300, 300),
        "margin": 20
    }
    igraph.plot(pes_gra, target=axes, **visual_style)

    fig.set_size_inches(16, 12)
    fig.savefig(file_name, dpi=200)
    plt.close(fig)


def _node_lst(conn, sccs_label=''):
    """ nx formatted nodes from a single conn
    """
    return [(sccs_label + node, {'color': 'green'}) for node in conn]


def _edge_lst(conn, sccs_label=''):
    """ nx formatted edge from a single conn
    """
    return [(
       *[sccs_label + node for node in conn], {'color': 'green'})]


def _well_nodes(nodes):
    """ list of nodes that are wells
    """
    return [node for node in nodes if '+' not in node]


def _prod_nodes(nodes):
    """ list of nodes that are bimol
    """
    return [node for node in nodes if '+' in node]


def _prod_edges(edges):
    """ list of edges that are to bimols
    """
    prod_edges = []
    for edge in edges:
        if '+' in edge[0] or '+' in edge[1]:
            prod_edges.append(edge)
    return prod_edges


def _well_edges(edges):
    """ list of edges that are between wells
    """
    well_edges = []
    for edge in edges:
        if '+' not in edge[0] and '+' not in edge[1]:
            well_edges.append(edge)
    return well_edges


def sccs_graph(conn_lst, sccs_idx=None):
    """ Make a networkx graph argument
    """
    G = nx.Graph()
    sccs_label = ''
    if sccs_idx is not None:
        sccs_label = '{:g}_{:g}!'.format(*sccs_idx)
    for conn in conn_lst:
        G.add_nodes_from(_node_lst(conn, sccs_label))
        G.add_edges_from(_edge_lst(conn, sccs_label))
    # Need to set the y-values
    # networkx.set_node_attributes(nxg, atom_symbols, 'symbol')
    return G


def _cluster_label_dct(nodes):
    """ simplify the label for the plot
    """
    label_dct = {}
    for node in nodes:
        label = node.split('!')[0].replace('_', ':')
        if (label not in label_dct.values() and
                '+' not in node.split('!')[1]):
            label_dct[node] = label
        else:
            label_dct[node] = ''
    # Loop again incase any sccs have no wells
    for node in nodes:
        label = node.split('!')[0].replace('_', ':')
        if label not in label_dct.values():
            label_dct[node] = label
    return label_dct


def _label_dct(nodes):
    """ simplify the label for the plot
    """
    label_dct = {}
    n = 14
    for node in nodes:
        label = node.split('!')
        if len(label) == 1:
            label = label[0]
        else:
            label = label[1]
        label_lst = label.split('+')
        split_label_lst = []
        for lab in label_lst:
            lab = '\n'.join([lab[i:i+n] for i in range(0, len(lab), n)])
            split_label_lst.append(lab)
        label_dct[node] = '\n+\n'.join(split_label_lst)

    return label_dct


def show_sccs(G, plt_title=None, save=False):
    """ show/save sccs graph
    """
    plt.figure(figsize=(10, 10), dpi=80)
    pos = nx.spring_layout(G, seed=3113794652)

    options = {
        "edgecolors": "black", "node_size": 3000,
        "alpha": 0.9, "linewidths": 3}
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=_well_nodes(G.nodes), node_color="tab:cyan", **options)

    options = {
        "edgecolors": "black", "node_size": 1800,
        "alpha": 0.6, "linewidths": 3}
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=_prod_nodes(G.nodes), node_color="tab:red", 
        node_shape="s", **options)

    nx.draw_networkx_edges(
        G, pos, edgelist=_well_edges(G.edges), width=4.0,
        alpha=0.4, edge_color='tab:cyan')
    nx.draw_networkx_edges(
        G, pos, edgelist=_prod_edges(G.edges), width=2.0,
        alpha=0.6, edge_color='tab:red')
    nx.draw_networkx_labels(
        G, pos, font_size=8, labels=_label_dct(G.nodes), font_color="black")

    if plt_title:
        plt.title(plt_title)
    if save:
        pass
    else:
        plt.show()


def _cross_conn(conns_a, conns_b, a_idx, b_idx):
    """ find the connections between bimol products and wells
        of two PESes
    """
    conn_lst = ()
    a_nodes = []
    b_nodes = []
    lab_a_nodes = []
    lab_b_nodes = []
    sccs_label_a = '{:g}_{:g}!'.format(*a_idx)
    sccs_label_b = '{:g}_{:g}!'.format(*b_idx)
    for conn in conns_a:
        a_nodes.extend(_node_lst(conn))
        lab_a_nodes.extend(_node_lst(conn, sccs_label_a))
    for conn in conns_b:
        b_nodes.extend(_node_lst(conn))
        lab_b_nodes.extend(_node_lst(conn, sccs_label_b))
    for i, a_node in enumerate(a_nodes):
        a_node = a_node[0].split('+')
        for j, b_node in enumerate(b_nodes):
            b_node = b_node[0].split('+')
            if len(a_node) > len(b_node):
                if b_node[0] in a_node:# or b_node[0].replace('A0', 'A1') in a_node or b_node[0].replace('B0', 'B1') in a_node:
                    conn_lst += (
                        (lab_a_nodes[i][0], lab_b_nodes[j][0]),)
            elif len(b_node) > len(a_node):
                if a_node[0] in b_node:# or a_node[0].replace('A0', 'A1') in b_node or a_node[0].replace('B0', 'B1') in b_node:
                    conn_lst += (
                        (lab_a_nodes[i][0], lab_b_nodes[j][0]),)
            elif len(a_node) == 1 and len(b_node) == 1:
                # if a_node[0] not in b_node and (a_node[0].replace('A0', 'A1') in b_node or a_node[0].replace('B0', 'B1') in b_node):
                #     conn_lst += (
                #         (lab_a_nodes[i][0], lab_b_nodes[j][0]),)
                #elif a_node == b_node:
                if a_node == b_node:
                    conn_lst += (
                        (lab_a_nodes[i][0], lab_b_nodes[j][0]),)
    return conn_lst


def show_pes(G_lst, g_conn_lst, idx_lst, save=False, connect_enants=True):
    """ show/save graph of mechanism
    """
    cross_ccs_conns = ()
    #chosen = ((0, 0), (1, 1), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0), (7, 0))
    chosen = ((0, 0), (1, 1), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0), (7, 1), (8, 0), (9, 0), (10, 0), (11, 0), (12, 0), (13, 0), (14, 0), (15, 0), (16, 0), (17, 0), (18, 0), (19, 0), (20, 0), (21, 0), (22, 0), (23, 0), (24, 0))
    print('idx list', idx_lst)
    for idx_a, (ia, ja) in enumerate(idx_lst):
        for idx_b, (ib, jb) in enumerate(idx_lst):
            #if ia > ib:
            if True:
                # if ia in [9, 10, 11, 12, 13, 14]:
                #if ia in [2, 3, 4]:
                #    continue
                # if ib in [9, 10, 11, 12, 13, 14]:
                #if ib in [2, 3, 4]:
                #    continue
                conns_a = g_conn_lst[idx_a]
                conns_b = g_conn_lst[idx_b]
                cross_conns = _cross_conn(
                    conns_a, conns_b, (ia, ja), (ib, jb))
                cross_ccs_conns += cross_conns
    full_G = nx.Graph()
    # for interactive html
    # new_G = Network(height="750px", width="100%", select_menu=True, cdn_resources='remote')
    for G in G_lst:
        full_G.update(G)

    for conn in cross_ccs_conns:
        full_G.add_edges_from(_edge_lst(conn))
        print('added cross connection', conn)
    colors = [
        #'tab:red', 'tab:blue', 'tab:orange', 'tab:pink',
        #'tab:green', 'tab:purple', 'tab:cyan', 'tab:olive',
        #'tab:gray', 'tab:brown', 'tomato', 'lightgreen',
        #'red', 'blue', 'orange', 'pink',
        #'green', 'slateblue', 'lightseagreen', 'mediumvioletred',
        #'dodgedblue', 'purple', 'tomato', 'lightgreen',
        'red', 'lightseagreen', 'lime', 'chartreuse',
        'green', 'crimson', 'dodgerblue', 'seagreen',
        'darkseagreen', 'mediumvioletred', 'steelblue', 'lightgreen',
        'slateblue', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black'
        ]
    edge_colors = [
        #'tab:red', 'tab:blue', 'tab:orange', 'tab:pink',
        #'tab:green', 'tab:purple', 'tab:cyan', 'tab:olive',
        #'tab:gray', 'tab:brown', 'tomato', 'lightgreen', 
        'red', 'blue', 'orange', 'pink',
        'green', 'purple', 'cyan', 'olive',
        'gray', 'brown', 'tomato', 'lightgreen', 
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black',
        'black', 'black', 'black', 'black', 'black', 'black'
        ]
    plt.figure(figsize=(14, 14), dpi=80)
    pos = nx.spring_layout(full_G, seed=3113794652)
    for idx, G in enumerate(G_lst):
        i, j = idx_lst[idx]
        #if (i, j) in chosen:
        if True:
            alpha = 0.9
        else:
            alpha = 0.2
        options = {
            "edgecolors": edge_colors[j], "node_size": 600, "alpha": alpha,
            "linewidths": 3}
        nx.draw_networkx_nodes(
            full_G, pos,
            nodelist=_well_nodes(G.nodes), node_color=colors[i], **options)

        options = {
            "edgecolors": edge_colors[j], "node_size": 200, "alpha": alpha,
            "linewidths": 3, 'node_shape':'s'}
        nx.draw_networkx_nodes(
            full_G, pos,
            nodelist=_prod_nodes(G.nodes), node_color=colors[i], **options)

        nx.draw_networkx_edges(
            full_G, pos, edgelist=_well_edges(G.edges), width=2.0, alpha=alpha,
            edge_color=colors[i])
        nx.draw_networkx_edges(
            full_G, pos, edgelist=_prod_edges(G.edges), width=1.0, alpha=alpha,
            edge_color=colors[i])
        # interactive html
        # options = {
        #     "edgecolors": edge_colors[j], "node_size": 600, "alpha": alpha,
        #     "linewidths": 3}
        # well_node_lst = _well_nodes(G.nodes)
        # new_G.add_nodes(
        #     well_node_lst,
        #     color=[colors[i]]*len(well_node_lst))
        # prod_node_lst = _prod_nodes(G.nodes)
        # new_G.add_nodes(
        #     prod_node_lst,
        #     shape=['square'] *len(prod_node_lst),
        #     color=[colors[i]]*len(prod_node_lst))
        # well_edge_lst = _well_edges(G.edges)
        # new_G.add_edges(
        #     well_edge_lst)
        #     #color=[colors[i]]*len(well_edge_lst))
        # prod_edge_lst = _prod_edges(G.edges)
        # new_G.add_edges(
        #     prod_edge_lst)
        #     #color=[colors[i]]*len(prod_edge_lst))
    chosen_edges = ()
    other_edges = ()
    for cross_conn in cross_ccs_conns:
        ia, ja = cross_conn[0].split('!')[0].split('_')
        ib, jb = cross_conn[1].split('!')[0].split('_')
        # if (int(ia), int(ja),) in chosen and (int(ib), int(jb),) in chosen:
        if True:
            chosen_edges += (cross_conn,)
        else:
            other_edges += (cross_conn,)
    nx.draw_networkx_edges(
        full_G, pos, edgelist=chosen_edges,
        width=2.0, alpha=.95, edge_color='black')
    nx.draw_networkx_edges(
        full_G, pos, edgelist=other_edges,
        width=1.0, alpha=.05, edge_color='black')
    # nx.draw_networkx_labels(full_G, pos, font_size=8, font_color="black")
    # nx.draw_networkx_labels(
    #     full_G, pos, font_size=8,
    #     labels=_label_dct(full_G.nodes), font_color="black")
    # nx.draw_networkx_labels(
    #     full_G, pos, font_size=28,
    #     labels=_cluster_label_dct(full_G.nodes), font_color="black")

    # interactive html
    # new_G.add_edges(
    #     chosen_edges)
    # new_G.add_edges(
    #     other_edges)

    if save:
        pass
    else:
        plt.show()
        # interactive html
        # new_G.show_buttons()
        # new_G.show('nx.html')


def pes_graph2(conn_lst, ene_dct=None, label_dct=None,
               file_name='surface.html'):
    """ Make a networkx graph argument
    """

    ene_dct, conn_lst, rgts_lst = _format(ene_dct, conn_lst, label_dct)

    nxg = networkx.Graph()
    nxg.add_nodes_from(rgts_lst)
    nxg.add_edges_from(conn_lst)
    # Need to set the y-values
    # networkx.set_node_attributes(nxg, atom_symbols, 'symbol')

    net = Network('3000px', '3000px')
    net.from_nx(nxg)
    # net.enable_physics(True)
    net.set_options("""
        "nodes": {
            "fixed": {
                "y": true
    }}
    """)
    net.show_buttons(filter_=['nodes', 'edges', 'physics'])
    net.save_graph(file_name)


def _format(ene_dct, conn_lst, label_dct):
    """ a
    """

    def _relabel(ene_dct, conn_lst, label_dct):
        """ relabel the graph
        """
        if ene_dct is not None:
            _ene_dct = {label_dct[name]: ene for name, ene in ene_dct.items()}
        else:
            _ene_dct = None

        _conn_lst = ()
        for conn in conn_lst:
            _conn_lst += ((label_dct[conn[0]], label_dct[conn[1]]),)

        return _ene_dct, _conn_lst

    def _rgts_lst(conn_lst):
        """ Get the full list of reagents
        """
        lst = ()
        for conn in conn_lst:
            lst += conn
        return automol.util.remove_duplicates_with_order(lst)

    if label_dct is not None:
        ene_dct, conn_lst = _relabel(ene_dct, conn_lst, label_dct)
    rgts_lst = _rgts_lst(conn_lst)

    return ene_dct, conn_lst, rgts_lst
