#Chase

import numpy as np
import pandas as pd
import pubchempy as pcp  # VPN needs to be off if using one
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.integrate import simpson
import networkx as nx

# https://www.nature.com/articles/s41597-020-00634-8#Sec3


def generate_chromo_df():
    """
    Generates the cleaned up chromophore database.
    Just load the saved pkl file instead.
    """
    chromo_db = pd.read_csv("chromo_db.csv")

    top_solvents = (
        chromo_db["Solvent"].sort_values(ascending=False).value_counts()[:10].index
    )

    solvent_names = {
        smiles: pcp.get_compounds(smiles, namespace="smiles")[0].iupac_name
        for smiles in top_solvents
    }

    chromo_db_filtered = chromo_db[chromo_db["Solvent"].isin(top_solvents)]

    cols_to_drop = [
        "Tag",
        "Lifetime (ns)",
        "Quantum yield",
        "log(e/mol-1 dm3 cm-1)",
        "abs FWHM (cm-1)",
        "emi FWHM (cm-1)",
        "abs FWHM (nm)",
        "emi FWHM (nm)",
        "Molecular weight (g mol-1)",
        "Reference",
    ]

    chromo_db_filtered.drop(columns=cols_to_drop, inplace=True)
    chromo_db_filtered = chromo_db_filtered.reset_index().drop(columns="index")
    chromo_db_filtered["Solvent_iupac"] = chromo_db_filtered["Solvent"].map(
        solvent_names
    )
    chromo_db_filtered.to_pickle("chromo_db_filtered.pkl")


def example_overlap_plot():

    chromo_db = pd.read_pickle("chromo_db_filtered.pkl")

    # get min and max from db
    x_min = (
        min(
            chromo_db["Absorption max (nm)"].min(), chromo_db["Emission max (nm)"].min()
        )
        - 100
    )
    x_max = (
        max(
            chromo_db["Absorption max (nm)"].max(), chromo_db["Emission max (nm)"].max()
        )
        + 100
    )
    x = np.linspace(x_min, x_max, 1000)

    # calculate normal distribution around max
    std_dev = 8
    chromo_1_abs = chromo_db["Absorption max (nm)"][0]
    chromo_2_em = chromo_db["Emission max (nm)"][0]

    chromo_1_dist = norm.pdf(x, loc=chromo_1_abs, scale=std_dev)
    chromo_2_dist = norm.pdf(x, loc=chromo_2_em, scale=std_dev)

    overlap = np.minimum(chromo_1_dist, chromo_2_dist)
    chromo_area = simpson(chromo_2_dist, x)
    overlap_area = simpson(overlap, x)
    percent_overlap = overlap_area / chromo_area

    plt.plot(x, chromo_1_dist, label="chromo 1 abs")
    plt.plot(x, chromo_2_dist, label="chromo 2 em")
    plt.title(f"Abs/Em Overlap ({percent_overlap:.2f}%)")
    plt.legend()


def generate_overlap_df(df, std_dev=10):

    chromo_series = df["Chromophore"]
    abs_series = df["Absorption max (nm)"]
    em_series = df["Emission max (nm)"]
    solv_series = df["Solvent_iupac"]

    # generate all possible pairs between chromophores and abs/em
    chromo_pairs = np.array(np.meshgrid(chromo_series, chromo_series)).T.reshape(-1, 2)
    abs_em_pairs = np.array(np.meshgrid(abs_series, em_series)).T.reshape(-1, 2)
    solvent_pairs = np.array(np.meshgrid(solv_series, solv_series)).T.reshape(-1, 2)

    # filter away pairs where chromo compared against itself
    # and keep where chromos are in the same solvent
    diff_chromo_same_solv = (chromo_pairs[:, 0] != chromo_pairs[:, 1]) & (
        solvent_pairs[:, 0] == solvent_pairs[:, 1]
    )
    chromo_pairs = chromo_pairs[diff_chromo_same_solv]
    abs_em_pairs = abs_em_pairs[diff_chromo_same_solv]
    solvent_pairs = solvent_pairs[diff_chromo_same_solv]

    # get min and max for x range
    x_min = min(abs_series.min(), em_series.min()) - 100
    x_max = max(abs_series.max(), em_series.max()) + 100
    x = np.linspace(x_min, x_max, 1000)

    # perform simpson integration of overlap between all pairs
    percent_overlap = []
    counter = 1
    for pair in abs_em_pairs:
        abs_dist = norm.pdf(x, loc=pair[0], scale=std_dev)
        em_dist = norm.pdf(x, loc=pair[1], scale=std_dev)
        overlap = np.minimum(abs_dist, em_dist)
        chromo_area = simpson(abs_dist, x)
        overlap_area = simpson(overlap, x)
        percent_overlap.append(overlap_area / chromo_area)
        print(f"pairs processed: {counter}/{len(abs_em_pairs)}", end="\r")
        counter += 1

    percent_overlap_df = pd.DataFrame()
    percent_overlap_df["Chromo A"] = chromo_pairs[:, 0]
    percent_overlap_df["Chromo B"] = chromo_pairs[:, 1]
    percent_overlap_df["A: abs (nm)"] = abs_em_pairs[:, 0]
    percent_overlap_df["B: em (nm)"] = abs_em_pairs[:, 1]
    percent_overlap_df["Solvent"] = solvent_pairs[:, 0]
    percent_overlap_df["Percent Overlap"] = percent_overlap

    percent_overlap_df.to_pickle("percent_overlap_df.pkl")


def create_chromo_graph(overlap_df, solvent, overlap_threshold=0.10):

    overlap_df = overlap_df["Percent Overlap"] >= overlap_threshold
    G = nx.from_pandas_edgelist(
        overlap_df,
        source=overlap_df["Chromo A"],
        target=overlap_df["Chromo B"],
        edge_attr=overlap_df["Percent Overlap"],
    )

    return G


if __name__ == "__main__":
    generate_chromo_df()
    chromo_db = pd.read_pickle("chromo_db_filtered.pkl")
    generate_overlap_df(chromo_db)
    overlap_db = pd.read_pickle("percent_overlap_df.pkl")
