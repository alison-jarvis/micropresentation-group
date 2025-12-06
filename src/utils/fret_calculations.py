import networkx as nx
from scipy import interpolate
import numpy as np
import math

# Constant
EPS_JREL = 1e-12

# Function to compute the spectral overlap integral from photochemCAD data
def compute_relative_J(
        lambda_em_nm,
        int_em,
        lambda_abs_nm,
        absorbance, grid_res=1.0):
    """
    Compute relative spectral overlap J_rel = integral F_D(lambda) * A_A(lambda) * lambda^4 dlambda
    Inputs: 
      - lambda_em_nm, int_em: donor emission wavelengths (nm) and intensities (arbitary units)
      - lambda_abs_nm, absorbance: acceptor wavelengths (nm) and absorbance A(lambda) (arbitrary units)
      - grid_res: Resolution of the wavelength grid for numerical integration, nominally 1 nm spacing over overlap.
    Returns:
      - J_rel (numeric, arbitrary units)
    """
    lam_min = max(min(lambda_em_nm), min(lambda_abs_nm))
    lam_max = min(max(lambda_em_nm), max(lambda_abs_nm))
    grid_nm = np.arange(lam_min, lam_max + grid_res, grid_res)

    f_em = interpolate.interp1d(lambda_em_nm, int_em, bounds_error=False, fill_value=0.0)
    f_abs = interpolate.interp1d(lambda_abs_nm, absorbance, bounds_error=False, fill_value=0.0)
    I_D = f_em(grid_nm)
    A_A = f_abs(grid_nm)

    delta = np.mean(np.diff(grid_nm))  # nm
    area = np.sum(I_D) * delta
    if area <= 0:
        raise ValueError("Donor emission has zero area after interpolation.")
    F_D = I_D / area

    # Use lambda in nm, include nm^4 factor - units are arbitrary so J_rel units also arbitrary
    J_rel = np.sum(F_D * A_A * (grid_nm**4)) * delta

    return J_rel

def generate_weight(C):
    def Jrel_to_weight(Jrel):
        e = 1 / (1 + (C/Jrel))
        return e
    return Jrel_to_weight

# Built the fret graph
def build_fret_graph(molecules, ems_abs):
    """
    Builds a directed graph where edge i to j represents fret efficiency of donor i exciting acceptor j
    Weight = -log(E(i to j))
    """

    # Directed graph
    G = nx.DiGraph()

    # Add nodes
    for m in molecules:
        G.add_node(m)

    Jrels = []

    # Add weighted edges
    for donor in molecules:
        for acceptor in molecules:
            if donor == acceptor:
                continue  # bc we don't want self FRET

            lam_em, inten_em = ems_abs[donor]['emission']['wavelength'], ems_abs[donor]['emission']['intensity']
            lam_abs, absorb = ems_abs[acceptor]['absorption']['wavelength'], ems_abs[acceptor]['absorption']['intensity']

            # Calculate J rel
            try:
                J_rel = compute_relative_J(
                    lam_em, inten_em,
                    lam_abs, absorb,
                    grid_res=1.0
                )
                Jrels.append(J_rel)
            except Exception as e:
                continue

    # For the arbitrary constant in E calculation
    J_rel_mean = np.nanmean(Jrels)

    # Specific function using the current constant
    Jrel_to_weight = generate_weight(J_rel_mean)

    # Add weighted edges
    for donor in molecules:
        for acceptor in molecules:
            if donor == acceptor:
                continue  # bc we don't want self FRET

            lam_em, inten_em = ems_abs[donor]['emission']['wavelength'], ems_abs[donor]['emission']['intensity']
            lam_abs, absorb = ems_abs[acceptor]['absorption']['wavelength'], ems_abs[acceptor]['absorption']['intensity']

            # Calculate J_rel
            J_rel = compute_relative_J(
                lam_em, inten_em,
                lam_abs, absorb,
                grid_res=1.0
            )

            if J_rel <= EPS_JREL or np.isnan(J_rel):
                continue

            # Convert to edge weight
            w = -math.log(Jrel_to_weight(J_rel))

            # Add edge
            G.add_edge(donor, acceptor, weight=w, J_rel=J_rel)

    return G
