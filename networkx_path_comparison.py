import pandas as pd
import numpy as np
from networkx.algorithms.shortest_paths import astar_path, dijkstra_path
from graph_helpers import create_chromo_graph, get_peak_distance_heuristic, print_path
from time import time

print("Loading Dataframes...")
try:
    chromo_df = pd.read_pickle("chromo_db_filtered.pkl")
except FileNotFoundError as e:
    raise FileNotFoundError("Please add the chromo_db_filtered.pkl file to this directory")

try:
    overlap_df = pd.read_pickle("percent_overlap_df.pkl")
except FileNotFoundError as e:
    raise FileNotFoundError("Please add the percent_overlap_df.pkl file to this directory")
overlap_df = overlap_df.dropna()


print("Constructing Graphs...")
# Construct a dictionary of FRET closeness graphs for each solvent
solvent_graphs = {}
for solvent in overlap_df['Solvent'].unique():
    solvent_graphs[solvent] = create_chromo_graph(overlap_df, solvent)

# Add -ln(overlap) as our edge weights for pathfinding
for graph in solvent_graphs.values():
    for edge in graph.edges(data=True):
        edge[2]["weight"] = -np.log(edge[2]["Percent Overlap"])

print("Done!")
print()
print('-'*30)
print('#'*30)
print('-'*30)
print()

"Solvent Graph Sizes: "
print('-'*30)
print({name: len(graph) for name, graph in solvent_graphs.items()})
print()
print('-'*30)
print('#'*30)
print('-'*30)
print()

# Test A* between two random donor/acceptors in chloroform
solvent = "chloroform"
start_donor = "O=C1NC(=O)c2ccc(N3CCCCC3)c3cccc1c23"
end_acceptor = "COc1ccc(-c2cc3c(s2)=CC2=C(C(F)(F)F)c4cc5sc(-c6ccc(OC)cc6)cc5n4[B-](F)(F)[N+]=32)cc1"

heuristic_function = get_peak_distance_heuristic(solvent, chromo_df)

start = time()
output_path_astar = astar_path(
    solvent_graphs[solvent],
    start_donor,
    end_acceptor,
    heuristic_function,
    weight='weight'
)
stop = time()
print("Path Found by A*: ")
print_path(output_path_astar, solvent_graphs[solvent])
print(f"runtime: {stop-start:.3f}")

print('-'*30)
print()

# Compare with dikstra's
start = time()
output_path_dijkstra = dijkstra_path(
    solvent_graphs[solvent],
    start_donor,
    end_acceptor,
    weight='weight'
)
stop = time()
print("Path Found by Dijkstra's: ")
print_path(output_path_dijkstra, solvent_graphs[solvent])
print(f"runtime: {stop-start:.3f}")
