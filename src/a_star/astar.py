import heapq
import networkx as nx
from typing import Callable, Hashable
from functools import lru_cache


"""Implementation based on this tutorial: https://www.geeksforgeeks.org/dsa/a-search-algorithm/,
    as well as wikipedia"""
def A_star(
        graph: nx.DiGraph,
        source: Hashable,
        target: Hashable,
        heuristic: Callable,
        cache_heuristic: bool=True
    ) -> list:
    """Uses the A* algorithm to find the shortest path between two nodes in a graph.


    Args:
        graph (nx.DiGraph): The graph to search
        source (Hashable): The node to begin the path
        target (Hashable): The node to end the path at
        heuristic (Callable): A function that takes in two nodes and returns a heuristic distance between them.
        cache_heuristic (bool, optional): Whether to cache calls to heuristic. Defaults to True.

    Returns:
        list: A list of nodes along the shortest path
    """    

    # Wrap heuristic in lru_cache if caching enabled
    if cache_heuristic:
        primary_heuristic = lru_cache(maxsize=None)(heuristic)
    else:
        primary_heuristic = heuristic

    # check if a valid smiles string in dataframe
    if not graph.has_node(source) or not graph.has_node(target):
        raise ValueError("Invalid source or target chromophore")

    # check if start and end are the same
    if source == target:
        raise ValueError("Source and target cannot be the same chromophore")

    # parent map for reconstructing path
    came_from = {}

    # dict to store the path details
    node_scores: dict[Hashable, float] = dict()

    node_scores[source] = 0.0 # g score
    came_from[source] = None # starting node comes from nowhere

    # initialize empty open list for nodes to be visited
    open_list = []

    # add current node to open list, this will contain the total cost f and the node id (smiles)
    heapq.heappush(open_list, (primary_heuristic(source, target), source))

    while open_list:
        # pop minimum value from heap
        _, current_id = heapq.heappop(open_list)

        # check if current node is the target:
        if current_id == target:
            path = []
            current = target
            while current is not None:  # source node has no parent
                path.append(current)
                current = came_from[current] # update current to parent
            return path[::-1]  # reverse list to get correct order

        for neighbor in graph.successors(current_id):
    
            # calculate graph distance to neighbor
            g_new = node_scores[current_id] + graph[current_id][neighbor]["weight"]

            if g_new < node_scores.get(neighbor, float('inf')):
                h_new = primary_heuristic(neighbor, target)
                f_new = g_new + h_new
                came_from[neighbor] = current_id
                node_scores[neighbor] = g_new

                # add to open list for further searching
                heapq.heappush(open_list, (f_new, neighbor))

    return None
