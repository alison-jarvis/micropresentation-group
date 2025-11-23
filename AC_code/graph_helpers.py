import networkx as nx
import pubchempy as pcp
import pandas as pd
import numpy as np
from typing import Callable
from functools import lru_cache

import math
import heapq


def create_chromo_graph(
    overlap_df: pd.DataFrame, solvent: str, overlap_threshold: float = 0.10
) -> nx.DiGraph:
    """Construct a Digraph for each solvent"""
    overlap_df = overlap_df[overlap_df["Percent Overlap"] >= overlap_threshold]
    overlap_df = overlap_df[overlap_df["Solvent"] == solvent]
    G = nx.from_pandas_edgelist(
        overlap_df,
        source="Chromo B",
        target="Chromo A",
        edge_attr="Percent Overlap",
        create_using=nx.DiGraph,
    )
    return G


def get_peak_distance_heuristic(solvent: str, chromo_df: pd.DataFrame) -> Callable:
    """
    Returns a function that, for a given solvent, calculates the distance
    in nm between emission an absorption peaks of donor and acceptor as a
    heuristic estimate of FRET yield. Note that this will not always lead
    to the shortest path if used as the herustic in A*, since the actual
    distance is not bounded below by it. In fact, it kinda sucks.
    """
    chromo_df_solvent = chromo_df.loc[chromo_df["Solvent_iupac"] == solvent]

    def heuristic(chromo_b, chromo_a):
        emission_wl = chromo_df_solvent.loc[
            chromo_df_solvent["Chromophore"] == chromo_b, "Emission max (nm)"
        ].iat[0]
        absorbtion_wl = chromo_df_solvent.loc[
            chromo_df_solvent["Chromophore"] == chromo_a, "Absorption max (nm)"
        ].iat[0]
        return abs(emission_wl - absorbtion_wl)

    return heuristic


### A bunch of stuff for printing out relatively pretty molecule strings in path traces ###
def abridge_smiles_str(smiles_str: str, max_len: int = 20) -> str:
    """Returns an abridge version of a smiles string."""
    if len(smiles_str) <= max_len:
        return smiles_str

    abridged_str = "...".join(
        [smiles_str[: max_len // 2], smiles_str[-(max_len - max_len // 2) :]]
    )
    return abridged_str


@lru_cache(maxsize=2048)
def smiles_to_uipac_name(smiles_str: str) -> str:
    uipac_name = pcp.get_compounds(smiles_str, namespace="smiles")[0].iupac_name
    if uipac_name:
        return abridge_smiles_str(uipac_name)
    else:
        return abridge_smiles_str(smiles_str)


def print_path(path: list[str], graph: nx.Graph) -> None:
    prev_node, cur_node = None, path[0]
    total_dist = 0
    for i, node in enumerate(path[1:]):
        prev_node = cur_node
        cur_node = node
        edge_weight = graph[prev_node][cur_node]["weight"]
        total_dist += edge_weight
        disp_prev, disp_cur = smiles_to_uipac_name(prev_node), smiles_to_uipac_name(
            cur_node
        )
        print(f"{i+1}. {disp_prev}-[{edge_weight:.2f}]->{disp_cur}")
    print(f"total distance: {total_dist:.2f}")
    print(f"effective yield: {np.exp(-total_dist):.2f}")


"""Implementation based on this tutorial: https://www.geeksforgeeks.org/dsa/a-search-algorithm/"""


class Chromo:
    def __init__(self):
        self.parent = 0
        self.f = float("inf")
        self.g = float("inf")
        self.h = 0


def A_star(graph: nx.Graph, source: str, target: str, heuristic) -> list:
    """Returns a list of nodes from shortest path"""

    # check if a valid smiles string in dataframe
    if not graph.has_node(source) or not graph.has_node(target):
        print("Invalid source or target chromophore")
        return

    # check if start and end are the same
    if source == target:
        print("Source and target cannot be the same chromophore")
        return

    # need a set for found nodes (chromophores)
    closed_list = set()

    # dict to store the path details
    chromo_details = {}

    # initialize the source node details, going to make a chromophore class
    source_node = Chromo()
    source_node.f = 0
    source_node.g = 0
    source_node.h = 0
    source_node.parent = None

    chromo_details[source] = source_node

    # initialize empty open list for nodes to be visited
    open_list = []

    # add current node to open list, this will contain the total cost f and the node id (smiles)
    heapq.heappush(open_list, (0.0, source))

    while open_list:
        # pop minimum value from heap
        _, current_id = heapq.heappop(open_list)

        # check if current node is the target:
        if current_id == target:
            path = []
            current = target
            while current is not None:  # source node has no parent
                path.append(current)
                current = chromo_details[current].parent  # update current to parent
            return path[::-1]  # reverse list to get correct order

        # check if current node id has already been visited
        if current_id in closed_list:
            continue

        closed_list.add(current_id)

        for neighbor in graph.successors(current_id):
            # only check neighbors to current node if not previously visited
            if neighbor in closed_list:
                continue

            # calculate costs for neighbor
            g_new = chromo_details[current_id].g + graph[current_id][neighbor]["weight"]
            h_new = heuristic(neighbor, target)
            f_new = g_new + h_new

            # if path to neighbor has lower cost then update costs
            if neighbor not in chromo_details or chromo_details[neighbor].f > f_new:

                # add to open list for further searching
                heapq.heappush(open_list, (f_new, neighbor))

                # add new node to chromo_details for path
                # if a lower cost path, value for key gets overwritten
                new_node = Chromo()
                new_node.f = f_new
                new_node.g = g_new
                new_node.h = h_new
                new_node.parent = current_id

                chromo_details[neighbor] = new_node

    print("No path found")
    return None
