import networkx as nx
import pubchempy as pcp
import pandas as pd
import numpy as np
from typing import Callable
from functools import lru_cache

def create_chromo_graph(overlap_df:pd.DataFrame, solvent:str, overlap_threshold:float=0.10) -> nx.DiGraph:
    '''Construct a Digraph for each solvent'''
    overlap_df = overlap_df[overlap_df["Percent Overlap"] >= overlap_threshold]
    overlap_df = overlap_df[overlap_df["Solvent"] == solvent]
    G = nx.from_pandas_edgelist(
        overlap_df,
        source="Chromo B",
        target="Chromo A",
        edge_attr="Percent Overlap",
        create_using=nx.DiGraph
    )
    return G

def get_peak_distance_heuristic(solvent:str, chromo_df:pd.DataFrame) -> Callable:
    '''
    Returns a function that, for a given solvent, calculates the distance
    in nm between emission an absorption peaks of donor and acceptor as a
    heuristic estimate of FRET yield. Note that this will not always lead
    to the shortest path if used as the herustic in A*, since the actual
    distance is not bounded below by it. In fact, it kinda sucks.
    '''
    chromo_df_solvent = chromo_df.loc[chromo_df["Solvent_iupac"]==solvent]
    def heuristic(chromo_b, chromo_a):
        emission_wl = chromo_df_solvent.loc[chromo_df_solvent["Chromophore"]==chromo_b, "Emission max (nm)"].iat[0]
        absorbtion_wl = chromo_df_solvent.loc[chromo_df_solvent["Chromophore"]==chromo_a, "Absorption max (nm)"].iat[0]
        return abs(emission_wl-absorbtion_wl)
    return heuristic

### A bunch of stuff for printing out relatively pretty molecule strings in path traces ###
def abridge_smiles_str(smiles_str: str, max_len: int=20) -> str:
    '''Returns an abridge version of a smiles string.'''
    if len(smiles_str) <= max_len:
        return smiles_str
    
    abridged_str = "...".join([smiles_str[:max_len//2], smiles_str[-(max_len-max_len//2):]])
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
        edge_weight = graph[prev_node][cur_node]['weight']
        total_dist += edge_weight
        disp_prev, disp_cur = smiles_to_uipac_name(prev_node), smiles_to_uipac_name(cur_node)
        print(f"{i+1}. {disp_prev}-[{edge_weight:.2f}]->{disp_cur}")
    print(f"total distance: {total_dist:.2f}")
    print(f"effective yield: {np.exp(-total_dist):.2f}")

