This repository contains an implementation of the algorithm described in:

 ### J. Simmons and T. Carrington, "Computing vibrational spectra using a new collocation method with a pruned basis and more points than functions: avoiding quadrature." The Journal of Chemical Physics 158, 144115 (2023) ###

which uses rectangular collocation with LU decomposition and a mapping procedure to construct and solve pruned tensor-product representations for multidimensional quantum-molecular vibrational problems.

The code is currently configured for a 12-D bilinearly coupled harmonic oscillator as a solvable reference problem, though is designed to be flexible for a range of potentials and collocation grids. The program uses BLAS, LAPACK, and ARPACK 96 to compute a selected subset of eigenvalues and eigenvectors. The Fortran source has been updated to use `iso_fortran_env` for explicit `real64` kind declarations; otherwise, the implementation retains the original Fortran structure of the algorithm.
