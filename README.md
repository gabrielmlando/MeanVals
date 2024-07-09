# MeanVals

This code is a companion to the paper [Computing Mean Values in the Deep Chaotic Regime, Phys. Rev. Lett. 132 (2024), 260401, by Gabriel M. Lando, Olivier Giraud and Denis Ullmo]. Please cite the original paper as

\begin{quote}
@article{PhysRevLett.132.260401,
  title = {Computing Quantum Mean Values in the Deep Chaotic Regime},
  author = {Lando, Gabriel M. and Giraud, Olivier and Ullmo, Denis},
  journal = {Phys. Rev. Lett.},
  volume = {132},
  issue = {26},
  pages = {260401},
  numpages = {6},
  year = {2024},
  month = {Jun},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevLett.132.260401},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.132.260401}
}
\end{quote}

## `QuantumKRS` module

This module computes quantum evolution and mean values for the kicked rotor system (KRS) using the methods described in Sec.~1 of the Supplemental Material. Evolution is computed in the angle representation, and we can switch to angular momentum representation using fast Fourier transforms. For wavefunctions the output actually consists of both angular momentum and angle grids used in the Fourier grid scheme, together with wavefunctions in the corresponding representations. Mean values are computed according to (S.12).

## `SemiclassicalMeanValueKRS` module

### `get_Lk` function

This function performs exactly what is described in Subsec.~IIIA of the Supplemental Material

### `sorters.jl` file

The main formula in `get_means_osc` assumes that all filaments are indexed according to the same rule. Thus, after the filaments are obtained by `get_Lk`, we choose the following ordering: All filaments are ordered from left to right, from top to bottom. That is, the first element of every filament will always be its leftmost intersection with the detector, and the first filament will always be above the others in the angular momentum axis. This is in line with Fig.~S1.

## `HermanKlukKRS` module

This module provides a straightforward calculation of Herman-Kluk wavefunctions for the KRS, and should be easy to read even for the reader not familiarized with the use of initial value representations. Wavefunctions are also outputed in angular momentum and angle representations, together with the corresponding grids (just as `QuantumKRS`). The module itself does not compute mean values; instead, we compute all wavefunctions in parallel (using a cluster) and save them to `.txt` files, each having a size of approximately `1Mb` for a grid with $2^{13}$ points. For computing mean-values, we use (S.12) and the same recipe employed in `QuantumKRS`.
