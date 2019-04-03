# Cluster
DGB Implementation for analyzing a single water molecule within a water cluster under the Local Monomer Approximation. 
Solves the generalized eigenvalue problem using a GDB basis, overlap and kinetic matricies computed analytically, potential computed numerically with QMC-methods. 
The dimensionality of the Lmon approximation is adjustable (Lmon1-9). 

## Generate
Example bash scripts for generating input files (to vary the Gaussian Width parameter) and submit to cluster via qsub. 

## src
DGB_cluster source code (Fortran 90).

## Data_Analysis
Example scripts for plotting eigenvalues as a function of alpha and eigenvalues as a function of sobol point.
Example plots are generated using xmgrace, make_data.sh will generate Nsobol and alpha plots. 

You generate as many eigenvalues as basis functions used.
For water ~3000 wavenumbers above the ground state corresponds to the vibrational modes of interest.
clean_alpha.sh and clean_eig.sh have parameters for setting the number of eigenvalues to plot. 

# Requirements
To run the code you need an xyz water configuration, a potential surface, and a quasi-random sequence generator. 

### Potentials
TIP4P and MBPOL are implemented currently, requires their source code to compile (TIP4P.f90, mbpol). 

### Sobol.f90 
Fortran-90 module for generating standard Sobol Sequences. 
http://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html

### XYZ Configuration
Potentials only intended for water calculations (can be modified if desired, check the cluster_module). 
cage.xyz ==> minimized water configuration for the cage-hexamer (TIP4P). 
