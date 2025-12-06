## __Spectral Pathfinding: Using A* to Navigate Chromophore Networks for FRET Design__

Ana Barrera-Jauregui: anabarrera_jauregui@berkeley.edu

Colman Bouton: clbouton@berkeley.edu

Alex Chase: achase95@berkeley.edu

Alison Jarvis: alison_jarvis@berkeley.edu

**To use this repo**

- Prepare the pickle spectral database and run the spectral-overlap code to:
    - visualize the donor/acceptor spectral pairs
    - compute the r_rel for all ordered pairs

- Build the chromophore graph from the overlaps and assign edge weights

- Apply A* to find the optimal energy-transfers and compare with the graphs; given the start (donor) and finish (acceptor) node.

- Analyse the choosen chromophore graph to predict the energy-transfer pathways


**A Pathfiding on FRET Spectral Overlaps**

This project combines Förster resonance energy transfer (FRET) spectroscopy with graph-based A* pathfinding. Nodes represent chromophores, in our case natural chlorophylls, edges represent possible donor–acceptor energy-transfer steps, and edge weights are derived from spectral overlap integrals between donor emission and acceptor absorption spectra.

In this project: 
- the graph is defined by the chromophores and overlap-based edge weights

- the goal is to find a path from an initial donor to a final acceptor that optimized some property.

Project Overview
Load chlorophyll (or other chromophore) spectra from a dataset.

- Load chlorophyll spectra from the dataset.

- For each ordered donor-acceptor pair, compute a relative spectral overlap *J_rel*

- Use these overlaps to define the edge weights in the graphs

- Run A* search to find the optimal paths through a chromophore network that can maximize the energy-transfer while optimizing the cost efficiveness when deriving said overlaps. 


***FRET basics***


Förster resonance energy transfer (FRET) describes non-radiative energy transfer from an excited donor chromophore to a nearby acceptor via dipole–dipole coupling. The efficiency of FRET depends on distance, relative orientation, donor quantum yield, and the spectral overlap between donor emission and acceptor absorption. FRET Efficiency:

$$\( E = 1 + \left( \frac{R}{r} \right)^6 \)$$

where r is the donor-acceptor distance and Ro is the Forster radius. 

***From overlap to graph weights***

To use FRET-like overlaps in the graph:

*NODES* = Chromophore

*DIRECTED EDGES* = In each ordered donor-acceptor pair (D,A), we had to compute Jrel(D -> A )

*EDGE COST EFFICIENCY* = inverse of the overlap 

                **Larger Overlap  ⇒ Smaller cost  ⇒ Preferred Edges**

***A` search on the chromophore graph***

A* is the best-first search algorithm that is able to find the least-cost path from the graph by combining the cost so far with a heurisitc estimate of the cost to the goal. It maintains a priority queue ordered by **f(n) = g(n) + h(n)

- g(n) is the cost from the start of the node

- h(n) or C is a heuristic lower bound on the remaining cost  (NEED TO VERIFY)

A simple heuristic is to use a lower bound on the remaining path cost 
- the minimum outgoing edge cost among remaining nodes
- if h(n) = 0 


#### Team Contribution Statement

Ana Barrera-Jauregui: Attempted excitation graph representation from photochemCAD database to compare with absorbance/emission from selected chlorophyll strain. Introduced theory on FRET diagram on slides and replicated conceptual project scope to README.

Colman Bouton:

Alex Chase:
Contributed to initial project planning for implementation of FRET with the A* search algorithm. Processed chromophore dataset to generate chromophore interaction graph based on percent overlap of peak absorbance/emission. Wrote the code for the A* search algorithm to identify highest yielding multi-step FRET pathways between a source and target chromophore within the graph. Provided explanation of the A* algorithm in contrast with Dijkstra's algorithm for the micropresentation slides.

Alison Jarvis: 
Compiled and contextualized equations for FRET, created code to calculate FRET spectral overlap from emission and absorption data, found photochemCAD database and determined how to parse data from it and sort by solvent. Created corresponding graph representation of chlorophylls in the Et20 solvent and preliminary graph visualization. Ran Djikstra and A-star on this graph and analyzed results. Provided explanation of graph representation and A-star applied to graph for presentation. Researched photosynthetic FRET pathways and real-world applications. 
