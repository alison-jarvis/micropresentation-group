import numpy as np
import pandas as pd
import networkx as nx
import os

import os 
import pickle

database_path = './Natural Chlorophylls/'
abs_files, ems_files = get_absorption_emission_files(database_path)

def load_spectrum(file_path):
    """
    Load a wavelength-dependent spectrum from a CSV file.
    Must have two columns: wavelength(nm), intensity.
    """

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Spectrum file not found: {file_path}")

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise ValueError(f"Could not read CSV file: {file_path}\n{e}")

    if df.shape[1] < 2:
        raise ValueError(f"Spectrum file must have 2 columns (nm, intensity): {file_path}")

    wavelength = df.iloc[:, 0].values
    intensity = df.iloc[:, 1].values

    return wavelength, intensity

def normalize_spectrum(intensity, wavelength):
    area = np.trapz(intensity, wavelength)
    if area == 0:
        raise ValueError("Spectrum has zero area and cannot be normalized.")
    return intensity / area


def compute_overlap(wavelength, emitter, absorber):
    """
    Generic overlap integral:
        J = ∫ emitter(λ) * absorber(λ) * λ^4 dλ
    """
    λ4 = wavelength ** 4
    integrand = emitter * absorber * λ4
    return np.trapz(integrand, wavelength)


def compute_excitation_overlap(wavelength, emission, excitation):
    """
    Optional overlap integral for: emission into  excitation
    whether the next chromophore can be excited by the previous emission.
    """
    integrand = emission * excitation
    return np.trapz(integrand, wavelength)


def construct_fret_graph(chromophores, use_excitation=False):
    """
    chromophores: dict of the form
    {
        'ChromophoreName': {
            'emission_file': 'path/to/em.csv',
            'absorption_file': 'path/to/abs.csv',
            OPTIONAL:
            'excitation_file': 'path/to/ex.csv'
        },
        ...
    }

    returns: networkx.DiGraph with edge weights = J_ij
    """

    G = nx.DiGraph()
    spectra = {}


    for name, files in chromophores.items():
        try: #error handing from files
            λ_em, F = load_spectrum(str(files["emission_file"]))
            λ_abs, eps = load_spectrum(str(files["absorption_file"]))
        except Exception as e:
            raise RuntimeError(f"Error loading spectra for {name}:\n{e}")

        # Interpolate absorption onto emission wavelength grid
        eps_interp = np.interp(λ_em, λ_abs, eps)
        F_norm = normalize_spectrum(F, λ_em)

        # optional excitation spectrum (can be remmoved)
        excitation = None
        if use_excitation:
            if "excitation_file" not in files:
                raise ValueError(f"Excitation spectrum required but missing for {name}.")
            λ_ex, ex = load_spectrum(files["excitation_file"])
            excitation = np.interp(λ_em, λ_ex, ex)

        spectra[name] = {
            "wavelength": λ_em,
            "emission": F_norm,
            "absorption": eps_interp,
            "excitation": excitation
        }

        G.add_node(name)

    # Compute all J_ij overlaps
    for donor in spectra:
        for acceptor in spectra:
            if donor == acceptor:
                continue

            data_D = spectra[donor]
            data_A = spectra[acceptor]

            # base FRET overlap
            J = compute_overlap(
                data_D["wavelength"],
                data_D["emission"],
                data_A["absorption"]
            )

            # (fixing) excitation-overlap mode, multiply terms
            if use_excitation:
                if data_A["excitation"] is None:
                    raise ValueError(f"Excitation missing for {acceptor}")
                J_ex = compute_excitation_overlap(
                    data_D["wavelength"],
                    data_D["emission"],
                    data_A["excitation"]
                )
                # Combine the two metrics
                J = J * J_ex  # could be weighted differently

            G.add_edge(donor, acceptor, weight=J)

    return G

# Pathfinding: Dijkstra and A*

def convert_to_cost_graph(G):
    #Convert weight=J to cost = -J for path *maximization*.
    #(ref AC/graph_helpers.py)
    H = nx.DiGraph()
    for u, v, d in G.edges(data=True):
        J = d["weight"]
        cost = -J  # maximize original J by minimizing -J
        H.add_edge(u, v, weight=cost)
    return H


def dijkstra_best_path(G, start, end):

    #Returns the os.path with np.maximum total spectral overlap (J).

    
    H = convert_to_cost_graph(G)
    try:
        path = nx.dijkstra_path(H, start, end, weight="weight")
    except nx.NetworkXNoPath:
        raise RuntimeError(f"No path found between {start} and {end}")

    return path


def astar_best_path(G, start, end):

    """
    (ref AC/graph_helpers.py)
    Custom A* heuristic:
    h(n) = 0 (safe but not informative)
    You can improve this using absorption/excitation patterns.
    """

    H = convert_to_cost_graph(G)

    def h(n):
        # Heuristic disabled for safety & correctness.
        return 0

    try:
        path = nx.astar_path(H, start, end, heuristic=h, weight="weight")
    except nx.NetworkXNoPath:
        raise RuntimeError(f"No path found between {start} and {end}")

    return path

if __name__ == "__main__":
    chromos = {
        "Molecule1": {
            "emission_file": "mol1_em.csv",
            "absorption_file": "mol1_abs.csv",
            "excitation_file": "mol1_ex.csv"
        },
        "Molecule2": {
            "emission_file": "mol2_em.csv",
            "absorption_file": "mol2_abs.csv",
            "excitation_file": "mol2_ex.csv"
        }
    }

    G = construct_fret_graph(chromos, use_excitation=True)

    print("\nSpectral Overlap Edges:")
    for u, v, d in G.edges(data=True):
        print(f"{u} into {v} : J = {d['weight']:.4e}")

    print("\nBest Dijkstra Path:", dijkstra_best_path(G, "Molecule1", "Molecule2"))
    print("Best A* Path:", astar_best_path(G, "Molecule1", "Molecule2"))
