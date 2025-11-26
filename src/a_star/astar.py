import heapq
import networkx as nx

class Chromo:
    def __init__(self):
        self.parent = 0
        self.f = float("inf")
        self.g = float("inf")
        self.h = 0


"""Implementation based on this tutorial: https://www.geeksforgeeks.org/dsa/a-search-algorithm/"""
def A_star(graph: nx.Graph, source: str, target: str, heuristic) -> list:
    """Returns a list of nodes from shortest path"""

    # check if a valid smiles string in dataframe
    if not graph.has_node(source) or not graph.has_node(target):
        raise ValueError("Invalid source or target chromophore")

    # check if start and end are the same
    if source == target:
        raise ValueError("Source and target cannot be the same chromophore")

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
