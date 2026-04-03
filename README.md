# Lieb-Liniger Bethe Ansatz

This repository has some basic scripts to compute the ground state energy and excitation spectrum of the Lieb-Liniger model using the Bethe ansatz formalism. It only has basic functionality including computing the ground state energies and densities in the following situations:

- Thermodynamic limit in a uniform system
- Finite system with periodic or hard wall boundaries
- Non-uniform system using a local density approximation

Additionally, the excitation spectrum can be computed for all of the above. A couple of examples with instructive comments can be found in `/scripts`.

## To-do
- check correctness of LDA and add unit tests
- add function to compute particle density for hard wall system