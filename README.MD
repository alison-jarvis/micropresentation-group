## __Spectral Pathfinding: Using A* to Navigate Chromophore Networks for FRET Design__

Ana Barrera-Jauregui: anabarrera_jauregui@berkeley.edu

Colman Bouton: clbouton@berkeley.edu

Alex Chase: achase95@berkeley.edu

Alison Jarvis: alison_jarvis@berkeley.edu


## Background

Förster resonance energy transfer (FRET) describes the energy transfer that can occur between a donor and acceptor chromophore. This transfer of energy depends on both the emission and absorbance profiles of both chromophores as well as their proximity. Common applications of FRET include determination of proteins conformations, protein-protein interactions, chemical sensors, and many more. FRET also occurs within nature, a notable example including chlorophyll and the energy transfers that occur to facilitate the production of energy. The PhotoChemCAD database provides spectral information for many different compounds relevant to photochemistry, including different chlorophylls. We utilized data from this site to produce a graph of different chlorophylls, and their interactions with one another via FRET. Nodes represent the individual chromophore variants, and edges represent the FRET efficiency between these molecules. High FRET efficiency (E) is desirable, so in order for this to equate to a low cost edge, the $-( \log (E))$ was taken. Additionally the graph is directed due to the required energy flow from a donor to an acceptor chromophore. The FRET efficiency is defined via the following equation. Where $r$ is the separation distance between the donor and the acceptor and $R_o$ is the Förster radius which describes the distance at which the energy transfer is 50% efficient.

$E = \frac{1}{1 + \left(\dfrac{r}{R_o}\right)^6}$

## Implementation of A* Search Algorithm

A* is an informed search, driven by heuristics, that is able to find the least-cost path from a source node to a target node. It maintains a priority queue with a cost function defined as $f(n) = g(n) + h(n)$ to guide the search. This algorithm is based upon Dijkstra's algorithm, which similarly utilizes a priority queue but where the cost function is defined as $f(n) = g(n)$. Due to the inclusion of a heuristic function $h(n)$ the average time complexity for A* is often better compared to Dijkstra's. Dijkstra's algorithm will check all possible nodes within the connected component of the source and target nodes. In contrast, A* uses a heuristic function to estimate the cost from the current node in the path to the target node. With a reasonably effective heuristic function, the algorithm will need to check far fewer child nodes for each parent node along the path compared to Dijkstra's. Although in the worst case, time complexity for both algorithms is $O(E\log E)$ and the space complexity is $O(V)$. However in practice, due to the guiding nature of the heuristic function, A* has an average effective time complexity of $O(b^d)$ where $b$ is a branching factor that corresponds with the number of child nodes explored per parent node, and $d$ is the depth of the minimum cost search path. This makes A* oftentimes much more efficient at solving for the lowest cost path

## A* for Finding Lowest Cost FRET Pathway between Chlorophyll Variants

In our case, the lowest cost FRET pathway is one that maximizes FRET efficiency between a source chlorophyll and target chlorophyll molecule. Within our provided example code, we show how A* is used to find the lowest cost, highest efficiency path from Chlorophyll F to Chlorophyll A. These are biologically significant because Chlorophyll F is known to result in maximum attenuation and Chlorophyll A activates the electron transport chain for energy production within plants. Therefore they are reasonable candidates for a source donor and a target acceptor. We found the highest efficiency multi-step pathway between these two chromophore included an intermediate energy transfer with Chlorophyll D. While narrow in scope, this example highlights the potential use for path finding algorithms like A* in regard to finding the most photochemical pathways in a variety of applications.

## Limitations and Next Steps

The biggest challenge we faced when implementing this algorithm was locating a suitable set of data to evalute it on. Specifically, to actually calculate the FRET efficiency between two chromophores we required the quantum yields of the donor. Then, the heuristic function defined as:

$$ -log(\frac{1}{1-\frac{C'}{Q_{d}J_{rel}}}) $$

would constitute an admissible, bounding estimate of the distance between the donor and target.

Additionally, these were also required to calculate an admissible heuristic function for the algorithm. While these do not exist for chlorophylls specificially, there are databases that provide access to the required properties for proteins. In the future, we could use these as a further proof of concept.


## Code Execution

To test the code, first install the packages specified in `requirements.txt`. Then, run the jupyter notebook found at `/src/natural_chlorophyll_example.ipynb`.


#### Team Contribution Statement

Ana Barrera-Jauregui: Attempted excitation graph representation from photochemCAD database to compare with absorbance/emission from selected chlorophyll strain. Introduced theory on FRET diagram on slides and replicated conceptual project scope to README.

Colman Bouton:
Derived/proposed edge weights of $-log(E_{fret})$. Developed utility script to assess $A*$ performance, accuracy against networkx implementation and Dijkstra's. Assessed the admissibility of proposed heuristic functions. Updated $A*$ implementation to cache heuristic function calls, revised $A*$ algorithm for readability, efficiency. Devised ultimate admissibility schema given quantum yield. Assisted with visualization of best path in chlorophyll graph.

Alex Chase:
Contributed to initial project planning for implementation of FRET with the A* search algorithm. Processed chromophore dataset to generate chromophore interaction graph based on percent overlap of peak absorbance/emission. Wrote the code for the A* search algorithm to identify highest yielding multi-step FRET pathways between a source and target chromophore within the graph. Provided explanation of the A* algorithm in contrast with Dijkstra's algorithm for the micropresentation slides. Provided background on FRET, A* search and biological applications for ReadMe.

Alison Jarvis: 
Compiled and contextualized equations for FRET, created code to calculate FRET spectral overlap from emission and absorption data, found photochemCAD database and determined how to parse data from it and sort by solvent. Created corresponding graph representation of chlorophylls in the Et20 solvent and preliminary graph visualization. Ran Djikstra and A-star on this graph and analyzed results. Provided explanation of graph representation and A-star applied to graph for presentation. Researched photosynthetic FRET pathways and real-world applications. 
