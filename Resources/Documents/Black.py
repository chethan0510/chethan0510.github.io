# Takes command-line arguments as input: run the program as python3 Black.py <input file>
# The input is a file containing the "adjacency list" of the form:
# [[<successors of node 1>], [<successors of node 2>],...,[<successors of node n-1>],[]]
# The output is a pdf file containing the animation (to the local folder), the (black) pebbling complexity and the length of
# the shortest sequence (both to the terminal).

import sys
import matplotlib.pyplot
import networkx
import math
import ast
from matplotlib.patches import FancyArrowPatch
from matplotlib.backends.backend_pdf import PdfPages

sys.setrecursionlimit(10000)  # increase recursion depth


# Returns the set of vertices that can be pebbled in a particular configuration (excluding vertices already pebbled)
def find_candidates(graph, configuration):
    candidates = []

    for vertex in graph.nodes():
        # Check if source or if all predecessors pebbled
        if not set(graph.predecessors(vertex)) - set(configuration):
            candidates.append(vertex)

    return list(set(candidates) - set(configuration))


# Generates all the valid pebbling sequences
def find_sequence(graph, limit, current_sequence, current_space):
    # current_space denotes the space-complexity of current_sequence
    if len(current_sequence) == 0:
        current_configuration = []
    else:
        current_configuration = current_sequence[len(current_sequence)-1]

    if list(graph.nodes())[len(graph.nodes())-1] in current_configuration:  # i.e., if the sink carries a pebble
        sequences.append(current_sequence)

    # Recurse through all strategies
    # Recursion 1: All possible cases of "remove pebble"
    for vertex in current_configuration:
        next_configuration = list(current_configuration)
        next_configuration.remove(vertex)
        next_configuration.sort()

        # Check if configuration already reached using fewer pebbles
        if next_configuration in configuration_table:
            index = configuration_table.index(next_configuration)
            previous_space = space_table[index]

            if previous_space > current_space:  # Current recursion takes less space
                space_table[index] = current_space  # Modify table

                next_sequence = list(current_sequence)  # Recurse
                next_sequence.append(next_configuration)
                find_sequence(graph, limit, next_sequence, current_space)

            else:  # Continue: configuration already reached using fewer pebbles
                continue
        else:  # Add and recurse
            configuration_table.append(next_configuration)  # Memoize
            space_table.append(current_space)

            next_sequence = list(current_sequence)  # Recurse
            next_sequence.append(next_configuration)
            find_sequence(graph, limit, next_sequence, current_space)

    # Recursion 2: All possible cases of "place pebble"
    if len(current_configuration) < limit:  # Recurse only if limit on number of pebbles not reached
        # Candidate for placing a pebble
        candidates = find_candidates(graph, current_configuration)

        for vertex in candidates:
            next_configuration = list(current_configuration)
            next_configuration.append(vertex)
            next_configuration.sort()

            # New
            # Check if configuration already reached using fewer pebbles
            if next_configuration in configuration_table:
                index = configuration_table.index(next_configuration)
                previous_space = space_table[index]
                if previous_space > max(current_space, len(next_configuration)):  # Current recursion takes less space
                    space_table[index] = max(current_space, len(next_configuration))  # Modify table

                    next_sequence = list(current_sequence)  # Recurse
                    next_sequence.append(next_configuration)
                    find_sequence(graph, limit, next_sequence, max(current_space, len(next_configuration)))

                else:  # Continue: configuration already reached using fewer pebbles
                    continue
            else:  # Add and recurse
                configuration_table.append(next_configuration)  # Memoize
                space_table.append(max(current_space, len(current_configuration)))

                next_sequence = list(current_sequence)  # Recurse
                next_sequence.append(next_configuration)
                find_sequence(graph, limit, next_sequence, max(current_space, len(next_configuration)))

            # OLD
            if next_configuration in configuration_table:
                index = configuration_table.index(next_configuration)
                previous_space = space_table[index]

                if previous_space > current_space + 1:
                    continue
            else:
                next_sequence = list(current_sequence)
                next_sequence.append(next_configuration)
                find_sequence(graph, limit, next_sequence, max(current_space, len(next_configuration)))


# Draws graph with curved edges
def draw_network(canvas, graph, position, configuration):
    for vertex in graph:  # Draw nodes
        if vertex in configuration:
            c = matplotlib.pyplot.Circle(position[vertex], radius=0.1, fill=True, label='1')
        else:
            c = matplotlib.pyplot.Circle(position[vertex], radius=0.1, fill=False, label='1')
        canvas.add_patch(c)  # Adds node to the plot
        # matplotlib.pyplot.text(position[vertex][0],position[vertex][1],vertex,fontsize=5) # To-do: number the nodes
        graph.node[vertex]['patch'] = c

    for (u, v) in graph.edges():
        n1 = graph.node[u]['patch']  # Co-ordinates of u
        n2 = graph.node[v]['patch']  # Co-ordinates of v

        if v == u+1:
            e = FancyArrowPatch(n1.center, n2.center, patchA=n1, patchB=n2, arrowstyle='->',
                                mutation_scale=10.0, alpha=0.5)
        else:
            e = FancyArrowPatch(n1.center, n2.center, patchA=n1, patchB=n2, arrowstyle='->',
                                connectionstyle='arc3,rad=0.25', mutation_scale = 10.0, alpha=0.5)

        canvas.add_patch(e)  # Adds edges to the plot

    return e


def find_adjacency_list(graph):
    adjacency_list = []

    for vertex in graph.nodes():
        adjacency_list.append(list(graph.successors(vertex)))

    return adjacency_list


# Main body
# Check for syntax
if len(sys.argv) < 2:
    print("Syntax: python3 Black.py <input file>")
    quit()

# Read adjacency list from the file on to a list
input_file = open(sys.argv[1], 'r')
input_list = ast.literal_eval(input_file.read())
input_file.close()

# Construct graph from adjacency list
input_graph = networkx.DiGraph()
for i in range(1, len(input_list)+1):
    input_graph.add_node(i)

for i in range(1, len(input_list)+1):
    for j in input_list[i-1]:
        input_graph.add_edge(i, j)

# Main: finds space complexity and the space-optimal pebbling for input_graph
# Run over all number of pebbles
for i in range(1, len(input_list)+1):
    # Memoize: all reached configurations stored in configuration_table, along with the pebbles used in space_table
    configuration_table = [[]]
    space_table = [0]
    sequences = []
    find_sequence(input_graph, i, [[]], 0)

    if len(sequences) > 0:
        shortest_sequence = sequences[0]
        for sequence in sequences:
            if len(sequence) < len(shortest_sequence):
                shortest_sequence = sequence

        print("Space complexity: ", i)
        print("Length of the sequence: ", len(shortest_sequence))
        print("Shortest sequence: ", shortest_sequence)
        break

# Dictionary containing the position of vertices
default_position = {}
for i in range(1, len(input_graph)+1):
    default_position[i] = [i, 0]

# Create .pdf for the animation
pdf_pages = PdfPages(sys.argv[1]+'.pdf')
fig = matplotlib.pyplot.figure(figsize=(8, 4), dpi=100)

# A page per configuration
for configuration in shortest_sequence:
    # Plot starts
    ax = matplotlib.pyplot.gca()
    draw_network(ax, input_graph, default_position, configuration)
    ax.autoscale()
    matplotlib.pyplot.axis('equal')
    matplotlib.pyplot.axis('off')

    # Done with the plot
    pdf_pages.savefig(fig)
    fig.clf()

# Write the PDF document to the disk
pdf_pages.close()
