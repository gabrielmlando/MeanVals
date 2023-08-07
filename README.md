# MeanVals

This code is a companion to the paper "Computing mean values in the deep chaotic regime", soon to appear in published form. 

## \texttt{QuantumKRS} module

This module computes quantum evolution and mean values for the kicked rotor system (KRS) using the methods described in Sec.~1 of the Supplemental Material. Evolution is computed in the angle representation, and we can switch to angular momentum representation using fast Fourier transforms. For wavefunctions the output actually consists of both angular momentum and angle grids used in the Fourier grid scheme, together with wavefunctions in the corresponding representations. Mean values are computed according to (S.12).

## \texttt{SemiclassicalMeanValueKRS} module

### \texttt{get\_Lk} function

This function performs exactly what is described in Subsec.~IIIA of the Supplemental Material

### \texttt{sorters.jl} file

The main formula in \texttt{get\_means\_osc} assumes that all filaments are indexed according to the same rule. Thus, after the filaments are obtained by \texttt{get\_Lk}, we choose the following ordering: All filaments are ordered from left to right, from top to bottom. That is, the first element of every filament will always be its leftmost intersection with the detector, and the first filament will always be above the others in the angular momentum axis. This is in line with Fig.~S1.

## \texttt{HermanKlukKRS} module

This module provides a straightforward calculation of Herman-Kluk wavefunctions for the KRS, and should be easy to read even for the reader not familiarized with the use of initial value representations. Wavefunctions are also outputed in angular momentum and angle representations, together with the corresponding grids (just as \texttt{QuantumKRS}). The module itself does not compute mean values; instead, we compute all wavefunctions in parallel (using a cluster) and save them to \texttt{.txt} files, each having a size of approximately \texttt{1Mb} for a grid with $2^{13}$ points. For computing mean-values, we use (S.12) and the same recipe employed in \texttt{QuantumKRS}.
