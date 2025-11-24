import numpy as np
import pandas as pd
import networkx as nx
#from scipy.integrate import simps

def load_spectrum(file_path):
    df = pd.read_csv(file_path)
    wavelength = df.iloc[:, 0].values
    intensity = df.iloc[:, 1].values
    return wavelength, intensity


def normalize_emission(emission):
    integral = simps(emission)
    return emission / integral


def spectral_overlap(wavelength, F_D, epsilon_A):
    lam4 = wavelength**4
    product = F_D * epsilon_A * lam4
    return simps(product, wavelength)


def construct_fret_graph(chromophores):
    """
    chromophores: dict of form
    {
        'Name': {
            'emission_file': "path/to/em.csv",
            'absorption_file': "path/to/abs.csv"
        },
        ...
    }
    networkx graph from fed edges through overlap
    """
    G = nx.DiGraph()

    spectra = {} #need
    for name, files in chromophores.items():
        λ_em, F = load_spectrum(files["emission_file"])
        λ_abs, eps = load_spectrum(files["absorption_file"])

        # Interpolate absorption to emission wavelengths
        eps_interp = np.interp(λ_em, λ_abs, eps)

        # Normalize emission
        F_norm = normalize_emission(F)

        spectra[name] = {
            "wavelength": λ_em,
            "emission": F_norm,
            "absorption": eps_interp
        }

        G.add_node(name)

    # D-A J_ij
    for donor in spectra:
        for acceptor in spectra:
            if donor == acceptor:
                continue

            data_D = spectra[donor]
            data_A = spectra[acceptor]

            J = spectral_overlap(
                data_D["wavelength"],
                data_D["emission"],
                data_A["absorption"]
            )

            G.add_edge(donor, acceptor, weight=J)

    return G

#testing 2 out of 3 times have gotten a timeout error, but need to debug
if __name__ == "__main__":
    chromos = {
        "Molecule1": {
            "emission_file": "mol1_emission.csv",
            "absorption_file": "mol1_absorption.csv"
        },
        "Molecule2": {
            "emission_file": "mol2_emission.csv",
            "absorption_file": "mol2_absorption.csv"
        },
        "Molecule3": {
            "emission_file": "mol3_emission.csv",
            "absorption_file": "mol3_absorption.csv"
        }
    }

    graph = construct_fret_graph(chromos)

    print("Directed edges with spectral overlap weights:")
    for u, v, d in graph.edges(data=True):
        print(f"{u} into {v} : J = {d['weight']:.4e}")
