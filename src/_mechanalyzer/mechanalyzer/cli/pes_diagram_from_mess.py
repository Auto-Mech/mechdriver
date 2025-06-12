import copy
import matplotlib.pyplot as plt
import networkx as nx

def parse_mess_file(file_path, remove_fake=True):
    """
    Parses a MESS input file to extract species names and their energies.

    :param file_path: Path to the MESS input file
    :type file_path: str
    :return: Dictionary with species names as keys and energies as values
    :rtype: dict
    """
    species_dict = {}
    connection_dict = {}
    current_species_name = None
    current_species_type = None
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            # Species names
            if line.startswith('Well ') or line.startswith('Barrier') or line.startswith('Bimolecular'):
                current_species_name = line.split()[1]  # Extract species name
                current_species_type = line.split()[0]
                if line.startswith('Barrier'):
                    if remove_fake:
                        reac = line.split()[2].replace('FakeB-', '').replace('FakeW-', '')
                        prod = line.split()[3].replace('FakeB-', '').replace('FakeW-', '')
                    else:
                        reac, prod = line.split()[2:4]
                    if reac != prod:
                        connection_dict[current_species_name] = (
                            reac, prod)
                    else:
                        current_species_type = 'Fake'
            if 'Fake' in line and not 'Barrier' in line and remove_fake:
                current_species_type = 'Fake'

            # ZeroEnergy for wells and barriers
            elif line.startswith('ZeroEnergy'):
                energy = float(line.split()[-1])  # Extract energy value
                if current_species_type in ['Well', 'Barrier']:
                    species_dict[current_species_name] = energy

            # GroundEnergy for bimolecular species
            elif line.startswith('GroundEnergy'):
                energy = float(line.split()[-1])  # Extract energy value
                if current_species_type == 'Bimolecular':
                    species_dict[current_species_name] = energy

    # Remove any barriers that connect points without energies (i.e., dummy)
    for barrier, stable_points in copy.deepcopy(connection_dict).items():
        for stable_point in stable_points:
            if stable_point not in species_dict:
                connection_dict.pop(barrier)
                species_dict.pop(barrier)
                print(f'Removing barrier {barrier}, likely b/c it involves a dummy')
                continue

    return species_dict, connection_dict

def nudge_nodes_iteration(graph, x_positions, y_positions, min_distance=1.0, max_distance=30.0, jumble=0):
    """
    Adjust x positions to minimize line crossings while considering both x and y spacing.
    The minimum spacing constraint applies only to non-neighboring species.

    :param graph: NetworkX graph representing the connections
    :type graph: nx.Graph
    :param x_positions: Initial x positions of the nodes
    :type x_positions: dict
    :param y_positions: y positions of the nodes
    :type y_positions: dict
    :param min_distance: Minimum Euclidean distance between non-neighboring nodes
    :type min_distance: float
    :return: Refined x positions
    :rtype: dict
    """
    refined_positions = x_positions.copy()
    adj_dct = {}
    new_adj_dct = {}
    sorted_positions = sorted(refined_positions.items(), key=lambda x: x[1])  # Sort by x position

    # check each species with the species closest to it on the x axis
    # (or jumbled positions away)
    for i in range(1, len(sorted_positions)):
        prev_species, prev_x = sorted_positions[i - 1]
        if i + jumble >= len(sorted_positions):
            jumble = 0
        curr_species, curr_x = sorted_positions[i + jumble]
        prev_y = y_positions[prev_species]
        curr_y = y_positions[curr_species]

        # if the species are neighbors check the max distance
        # else check the min distance
        if graph.has_edge(prev_species, curr_species):
            pass
            # distance = abs(curr_x - prev_x)
            # if abs(curr_x) > 4:
            #     direction = -1 if (curr_x) > 0 else 1
            #     refined_positions[curr_species] = curr_x + direction
            #     print(f"Adjusting {curr_species} position from {curr_x} to {refined_positions[curr_species]} to bring it closer to 0")
            # elif distance > max_distance:
            #     direction = -1 if (curr_x - prev_x) > 0 else 1
            #     refined_positions[curr_species] = curr_x + direction * abs((curr_x - prev_x)) / 3
            #     print(f"Adjusting {curr_species} position from {curr_x} to {refined_positions[curr_species]} to bring it closer to {prev_species}")
        else:
            distance = ((curr_x - prev_x) ** 2 + (curr_y - prev_y) ** 2) ** 0.5
            if distance < min_distance:
                direction = 1 if (curr_x - prev_x) > 0 else -1
                direction /= 4
                refined_positions[curr_species] = curr_x + direction * (min_distance ** 2 - (curr_y - prev_y) ** 2) ** 0.5
                adj_dct[curr_species] = direction * (min_distance ** 2 - (curr_y - prev_y) ** 2) ** 0.5

    # for any node that was adjusted, adjust any neighbor it has in the nudged
    # direction by half the amound
    for edges in graph.edges:
        if edges[0] in adj_dct:
            if refined_positions[edges[0]] < refined_positions[edges[1]] and adj_dct[edges[0]] > 0:
                refined_positions[edges[1]] += adj_dct[edges[0]] / 2
                new_adj_dct[edges[1]] = adj_dct[edges[0]] / 2
            elif refined_positions[edges[0]] > refined_positions[edges[1]] and adj_dct[edges[0]] < 0:
                refined_positions[edges[1]] += adj_dct[edges[0]] / 2
                new_adj_dct[edges[1]] = adj_dct[edges[0]] / 2
        elif edges[1] in adj_dct:
            if refined_positions[edges[1]] < refined_positions[edges[0]] and adj_dct[edges[1]] > 0:
                refined_positions[edges[0]] += adj_dct[edges[1]] / 2
                new_adj_dct[edges[0]] = adj_dct[edges[1]] / 2
            elif refined_positions[edges[1]] > refined_positions[edges[0]] and adj_dct[edges[1]] < 0:
                refined_positions[edges[0]] += adj_dct[edges[1]] / 2
                new_adj_dct[edges[0]] = adj_dct[edges[1]] / 2

    # another iteration for the next neighbors
    for edges in graph.edges:
        if edges[0] in new_adj_dct and edges[1] not in adj_dct:
            if refined_positions[edges[0]] < refined_positions[edges[1]] and new_adj_dct[edges[0]] > 0:
                refined_positions[edges[1]] += new_adj_dct[edges[0]] / 2
            elif refined_positions[edges[0]] > refined_positions[edges[1]] and new_adj_dct[edges[0]] < 0:
                refined_positions[edges[1]] += new_adj_dct[edges[0]] / 2
        elif edges[1] in new_adj_dct and edges[0] not in adj_dct:
            if refined_positions[edges[1]] < refined_positions[edges[0]] and new_adj_dct[edges[1]] > 0:
                refined_positions[edges[0]] += new_adj_dct[edges[1]] / 2
            elif refined_positions[edges[1]] > refined_positions[edges[0]] and new_adj_dct[edges[1]] < 0:
                refined_positions[edges[0]] += new_adj_dct[edges[1]] / 2

    return refined_positions


def initiate_graph(species_dict, connection_dict):
    """
    Initiates a graph based on species and their connections.

    :param species_dict: Dictionary with species names as keys and energies as values
    :type species_dict: dict
    :param connection_dict: Dictionary with species names as keys and their connections as values
    :type connection_dict: dict
    :return: Graph object representing the connections
    :rtype: nx.Graph
    :return: Dictionary with species names as keys and their degrees as values
    :rtype: dict
    """
    graph = nx.Graph()
    degrees = {}
    for species, connections in connection_dict.items():
        for connected_species in connections:
            if connected_species not in degrees:
                degrees[connected_species] = 1
            else:
                degrees[connected_species] += 1
            degrees[species] = 1
    for species, connections in connection_dict.items():
        for connected_species in connections:
            if '+' in connected_species and degrees[connected_species] < 2:
                weight = 5
            else:
                weight = 3
            graph.add_edge(species, connected_species, weight=weight)

    return graph, degrees


def set_color_palette(colors_on, wells=[]):
    """
    Sets the color palette for the plot.

    :param colors_on: Whether to use colors for wells and their connections
    :type colors_on: bool
    :return: List of colors
    :rtype: list
    """
    if colors_on:
        return [
            'blue', 'green', 'red', 'orange',
            'purple', 'yellow', 'pink',
            'cyan', 'brown', 'gray']
    else:
        return ['black'] * len(wells)


def determine_wells(degrees, well_threshold=2):
    """
    Determines which species to use as centered wells based on the degree of each species.

    :param degrees: Dictionary with species names as keys and their degrees as values
    :type degrees: dict
    :param well_threshold: Minimum number of connections for a species to be considered a well
    :type well_threshold: int
    :return: List of wells
    :rtype: list
    """
    wells = []
    for (species, degree) in sorted(degrees.items(), key=lambda x: x[1], reverse=True):
        if degree > well_threshold:
            wells.append(species)
        else:
            break
    return wells

def set_initial_positions_and_colors(graph, wells, degrees, colors_on=True):
    """
    Sets the initial positions and colors for the wells and species.

    :param graph: Graph object representing the connections
    :type graph: nx.Graph
    :param wells: List of wells
    :type wells: list
    :param degrees: Dictionary with species names as keys and their degrees as values
    :type degrees: dict
    :param colors: List of colors
    :type colors: list
    :param colors_on: Whether to use colors for wells and their connections
    :type colors_on: bool
    :return: Dictionary with species names as keys and their positions as values
    :rtype: dict
    :return: Dictionary with species names as keys and their colors as values
    :rtype: dict
    """
    colors = set_color_palette(colors_on, wells)
    pos = {}
    node_colors = {}
    offset = 0.3 # making no well at 0 helps some of the other steps later

    for i, well in enumerate(wells):
        node_colors[well] = colors[i]

    # this loops through the species, starting with the most connected
    # and then assigns them a position further from center from their closestly connected well
    # if its equally close to two wells, it will be placed in between them
    # and a color based on their closestly connected well
    for i, (species, _) in enumerate(
            sorted(degrees.items(), key=lambda x: x[1], reverse=True)):
        if species in wells:
            pos[species] = [(i - len(wells)/2 + offset)*2]
        else:
            lengths = nx.single_source_shortest_path_length(graph, species)
            closest_targets = [
                (node, dist) for node, dist in lengths.items() if node in wells]
            min_dist = min(dist for _, dist in closest_targets)
            closest_targets = [
                node for node, dist in closest_targets if dist == min_dist]
            if len(closest_targets) > 1:
                if not '+' in species:
                    pos[species]  = [(pos[closest_targets[0]][0] + pos[closest_targets[1]][0])/2]
                    node_colors[species] = 'black'
                else:
                    pos[species]  = [-5] if pos[closest_targets[0]][0] + pos[closest_targets[1]][0] < 0 else [5]
                    node_colors[species] = node_colors.get(closest_targets[0], 'gray')
            else:
                pos[species]  = [pos[closest_targets[0]][0] + pos[closest_targets[0]][0] * min_dist/5.]
                node_colors[species] = node_colors.get(closest_targets[0], 'gray')

    return pos, node_colors


def nudge_nodes_iteratively(
        graph, degrees, pos, species_dict, connection_dict, wells, scaling_factor=1.0,
        min_distance=1.0, max_distance=30.0, nudge_iterations=10):
    """
    Adjust x positions to minimize overlap while considering both x and y spacing.
    The minimum spacing constraint applies only to non-neighboring species.
    :param graph: NetworkX graph representing the connections
    :type graph: nx.Graph
    :param degrees: Dictionary with species names as keys and their degrees as values
    :type degrees: dict:
    :param pos: Initial x positions of the nodes
    :type pos: dict
    :param species_dict: Dictionary with species names as keys and energies as values
    :type species_dict: dict
    :param wells: List of wells
    :type wells: list
    :param min_distance: Minimum Euclidean distance between non-neighboring nodes
    :type min_distance: float
    :param max_distance: Maximum Euclidean distance between neighboring nodes
    :type max_distance: float
    :param nudge_iterations: Number of iterations for nudging
    :type nudge_iterations: int
    :param jumble: Jumble factor for adjusting positions
    :type jumble: int
    :return: Refined x positions
    :rtype: dict
    """
    # when there are two paths to the same product, force the product
    # to be placed at the end of the furthest out path
    # instead of in the middle of the two paths, where networkx puts it
    for i, (species, _) in enumerate(
            sorted(degrees.items(), key=lambda x: x[1], reverse=True)):
        if not '+' in species:
            continue
        connected_species = [key for key, value in connection_dict.items() if species in value]
        if len(connected_species) > 1:
            dominant = (
                pos[connected_species[0]][0]
                if abs(pos[connected_species[0]][0]) >  abs(pos[connected_species[1]][0])
                else pos[connected_species[1]][0])
            pos[species]  = [dominant + .5] if dominant > 0 else [dominant - .5]

    # if there is a bunch of interconnections between wells in the middle of the PES
    # make sure that any product paths are pushed out further than these interconnections
    max_well_x = max(val[0] for species, val in pos.items() if species in wells)
    min_well_x = min(val[0] for species, val in pos.items() if species in wells)
    for species, value in pos.items():
        if not '+' in species:
            continue
        if value[0] > 0 and value[0] < max_well_x:
            pos[species] = [2*value[0] + max_well_x]
            for ts, species_lst in connection_dict.items():
                if species in species_lst:
                    pos[ts] = [pos[ts][0] + max_well_x/2]
        elif value[0] < 0 and value[0] > min_well_x:
            pos[species] = [2*value[0] + min_well_x]
            for ts, species_lst in connection_dict.items():
                if species in species_lst:
                    pos[ts] = [pos[ts][0] + min_well_x/2]

    # now iteratively make small nudges if the distance between two species is too large
    # or too small, and then nudge the neighbors of those species
    x_positions = {
        species: pos[species][0] * scaling_factor
        for species in species_dict.keys()}
    plot_range = max(val[0] for val in pos.values()) - min(val[0] for val in pos.values())
    for i in range(nudge_iterations):
        x_positions = nudge_nodes_iteration(
            graph, x_positions, species_dict,
            min_distance=min_distance, max_distance=max_distance,
            jumble=i%4)
    return x_positions


def generate_plot(
        species_dict, connection_dict, wells, x_positions, node_colors,
        output_file="pes_diagram", format="svg", aspect_ratio=1, labels=True):
    """
    Generates a plot of the PES diagram.

    :param species_dict: Dictionary with species names as keys and energies as values
    :type species_dict: dict
    :param connection_dict: Dictionary with species names as keys and their connections as values
    :type connection_dict: dict
    :param wells: List of wells
    :type wells: list
    :param x_positions: Dictionary with species names as keys and their positions as values
    :type x_positions: dict
    :param node_colors: Dictionary with species names as keys and their colors as values
    :type node_colors: dict
    :param output_file: Name of the output figure file
    :type output_file: str
    :param format: Format of the output figure file (e.g., svg, png)
    :type format: str
    :param aspect_ratio: Aspect ratio of the output figure (width / height)
    :type aspect_ratio: float
    :param labels: Whether to label the species in the PES
    :type labels: bool
    """
    fig, axes = plt.subplots(1, 1, figsize=(8*aspect_ratio, 8))
    average_y = sum(species_dict.values()) / len(species_dict)
    y_range = max(species_dict.values()) - min(species_dict.values())

    # plot nodes
    for species, energy in species_dict.items():
        x = x_positions[species]
        y = energy
        if y > average_y:
            label_shift = y_range / 70
        else:
            label_shift = -y_range / 35
        axes.barh(
            y, width=0.4, left=x - 0.2,
            color=node_colors.get(species, 'gray'),
            edgecolor=node_colors.get(species, 'gray'),
            height=1, zorder=2)
        if labels:
            axes.text(x, y + label_shift, species, fontsize=10, ha='center', zorder=3)

    # plot edges
    for species, connections in connection_dict.items():
        x1, y1 = x_positions[species], species_dict[species]
        color=(
            'black' if all(
                connected_species in wells for connected_species in connections)
            else node_colors.get(species, 'gray'))
        for connected_species in connections:
            x2, y2 = x_positions[connected_species], species_dict[connected_species]
            if x1 < x2:
                x1adj = x1 + .2
                x2 -= .2
            else:
                x2 += .2
                x1adj = x1 - .2
            axes.plot([x1adj, x2], [y1, y2], color=color, alpha=.4, linewidth=2, zorder=1)
    # Step 5: Customize the plot
    axes.set_ylabel("Energy (kcal/mol)")
    axes.set_xticks([])
    plt.title("Potential Energy Surface (PES) Diagram")
    plt.tight_layout()
    plt.savefig(output_file + "." + format, format=format, dpi=300, bbox_inches='tight')

def main(
       input_file, well_threshold=2, colors_on=True,
        gravity=1, spring_iterations=20000, nudge_iterations=10,
        min_distance=1.0, max_distance=30.0, output_file="pes_diagram", format="svg",
        aspect_ratio=1, labels=True, remove_fake=True, shift_energy=True):
    """
    Plots a PES diagram based on species energies and connections.

    :param species_dict: Dictionary with species names as keys and energies as values
    :type species_dict: dict
    :param connection_dict: Dictionary with species names as keys and their connections as values
    :type connection_dict: dict
    :param well_threshold: Minimum number of connections for a species to be considered a well
    :type well_threshold: int
    :param colors_on: Whether to use colors for wells and their connections
    :type colors_on: bool
    :param gravity: Attraction factor for the spring layout
    :type gravity: int
    :param spring_iterations: Number of iterations for the spring layout algorithm
    :type spring_iterations: int
    :param nudge_iterations: Number of iterations for nudging the species to minimize overlap
    :type nudge_iterations: int
    :param min_distance: Minimum distance between non-neighboring nodes
    :type min_distance: float
    :param max_distance: Maximum distance between neighboring nodes
    :type max_distance: float
    :param output_file: Name of the output figure file
    :type output_file: str
    :param format: Format of the output figure file (e.g., svg, png)
    :type format: str
    :param aspect_ratio: Aspect ratio of the output figure (width / height)
    :type aspect_ratio: float
    :param labels: Whether to label the species in the PES
    :type labels: bool
    :param shift_energy: Whether to shift energies relative to the lowest well
    :type shift_energy: bool
    """
    species_dict, connection_dict = parse_mess_file(input_file, remove_fake=remove_fake)
    graph, degrees = initiate_graph(species_dict, connection_dict)
    wells = determine_wells(degrees, well_threshold)
    if shift_energy:
        min_well_energy = min(
            species_dict[species] for species in wells) if wells else 0
        species_dict = {
            species: energy - min_well_energy for species, energy in species_dict.items()}
    # determine initial x positions for nodes
    pos, node_colors = set_initial_positions_and_colors(
        graph, wells, degrees, colors_on=colors_on)
    plot_range = max(val[0] for val in pos.values()) - min(val[0] for val in pos.values())
    scaling_factor = 1.0

    # let spring layout determine the semifinal positions
    pos = nx.spring_layout(
        graph, pos=pos,
        k=plot_range/gravity/8, scale=plot_range/2,
        iterations=spring_iterations)
    x_positions = {
        species: pos[species][0] * scaling_factor
        for species in species_dict.keys()}

    # nudge nodes to minimize overlap and prevent weird
    # product placement
    if nudge_iterations > 0:
        x_positions = nudge_nodes_iteratively(
            graph, degrees, pos, species_dict, connection_dict,
            wells, scaling_factor=scaling_factor,
            min_distance=min_distance, max_distance=max_distance,
            nudge_iterations=nudge_iterations)

    generate_plot(
        species_dict, connection_dict, wells, x_positions, node_colors,
        output_file=output_file, format=format, aspect_ratio=aspect_ratio, labels=labels)

if __name__ ==  "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Parse MESS input files and extract species energies.")
    parser.add_argument(
        "--input_file",
        "-i",
        type=str,
        help="Path to the MESS input file",
        default="mess.inp")
    parser.add_argument(
        "--well_threshold",
        "-w",
        type=int,
        help="How many connections a species must have to be centered as well",
        default=2)
    parser.add_argument(
        "--colors_on",
        "-c",
        type=bool,
        help="True/False colorful PES, automatically makes each well and their connections a unique color",
        default=True)
    parser.add_argument(
        "--gravity",
        "-g",
        type=int,
        help="How much the species are pulled together in the spring layout",
        default=1)
    parser.add_argument(
        "--spring_iterations",
        "-s",
        type=int,
        help="How many iterations to run the spring layout algorithm",
        default=20000)
    parser.add_argument(
        "--nudge_iterations",
        "-n",
        type=int,
        help="How many iterations to nudge the species to minimize overlap",
        default=10)
    parser.add_argument(
        "--min_distance",
        "-d",
        type=float,
        help="Minimum distance between non-neighboring nodes, aka whats considered overlap",
        default=1.0)
    parser.add_argument(
        "--labels",
        "-l",
        type=bool,
        help="True/False whether to label the species in the PES",
        default=True)
    parser.add_argument(
        "--output_file",
        "-o",
        type=str,
        help="Name of the output figure file",
        default="pes_diagram")
    parser.add_argument(
        "--format",
        "-f",
        type=str,
        help="Format of the output figure file (e.g., svg, png)",
        default="svg")
    parser.add_argument(
        "--aspect_ratio",
        "-a",
        type=float,
        help="Aspect ratio of the output figure (width / height)",
        default=1.0)
    parser.add_argument(
        "--shift_energy",
        "-e",
        type=bool,
        help="Whether to shift energies relative to the lowest well",
        default=True)
    args = parser.parse_args()

    main(
        input_file=args.input_file,
        well_threshold=args.well_threshold,
        colors_on=args.colors_on,
        gravity=args.gravity,
        spring_iterations=args.spring_iterations,
        nudge_iterations=args.nudge_iterations,
        min_distance=args.min_distance,
        max_distance=30.0,
        output_file=args.output_file,
        format=args.format,
        aspect_ratio=args.aspect_ratio,
        labels=args.labels,
        shift_energy=args.shift_energy)
