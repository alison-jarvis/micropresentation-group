** A Pathfiding on FRET Spectral Overlaps**

This project combines Förster resonance energy transfer (FRET) spectroscopy with graph-based A* pathfinding. Nodes represent chromophores, in our case natural chlorophylls, edges represent possible donor–acceptor energy-transfer steps, and edge weights are derived from spectral overlap integrals between donor emission and acceptor absorption spectra.


Project Overview
Load chlorophyll (or other chromophore) spectra from a dataset.

- Load chlorophyll spectra from the dataset.

- For each ordered donor-acceptor pair, compute a relative spectral overlap *J_rel*

- Use these overlaps to define the edge weights in the graphs

- Run A* search to find the optimal paths through a chromophore network that can maximize the energy-transfer while optimizing the cost efficiveness when deriving said overlaps. 


***FRET basics***


Förster resonance energy transfer (FRET) describes non-radiative energy transfer from an excited donor chromophore to a nearby acceptor via dipole–dipole coupling. The efficiency of FRET depends on distance, relative orientation, donor quantum yield, and the spectral overlap between donor emission and acceptor absorption. FRET Efficiency:

\( E = 1 + \left( \frac{R}{r} \right)^6 \)


where r is the donor-acceptor distance and Ro is the Forster radius. 
