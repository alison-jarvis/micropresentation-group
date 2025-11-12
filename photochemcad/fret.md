## Spectral Overlap Integral

The overlap of donor emission with acceptor absorption is defined by the spectral overlap integral, $J$. 

$$
J = \int_0^\infty F_D (\lambda) \epsilon _A (\lambda) \lambda^4\, d\lambda
$$

Where we have the following terms:

- $F_D$ : normalized donor emission spectrum (function of $\lambda$), unitless
- $\epsilon _A$ : acceptor molar extinction coefficient (function of $\lambda$), units of $M^{-1} cm^{-1}$
- $\lambda$ : wavelength, units of $cm$

We don't have the molar extinction coefficient. From photochemCAD we get the absorbance optical density. The calculation between these two follows the Beer-Lambert law:

$$A(\lambda) = \epsilon (\lambda) c l$$

Where $A(\lambda)$ is the absorbance optical density, $c$ is concentration, and $l$ is path length. Not sure where to find these values other than maybe the papers. If we assume they are consistent / constant, however, we can also assume that:

$$A(\lambda) \propto \epsilon (\lambda)$$

So we can calculate not an absolute absorbance but a relative one, since if $A(\lambda) \propto \epsilon (\lambda)$, it also follows that $J \propto J_{rel}$ where $J_{rel}$ is:

$$
J_{rel} = \int_0^\infty F_D (\lambda) A (\lambda) \lambda^4\, d\lambda
$$

Therefore, we can calculate an approximation of the spectral overlap using the formula for $J_{rel}$ and numerical integration of the donor emission ($F_D (\lambda)$) and the acceptor absorption ($A (\lambda)$). 