# openMPparallelExercise
Serial code for openMP exercise for parallel programming course

The absorption spectra of disordered linear molecular aggregates can be described with a so-called Frenkel Hamiltonian. In its simplest form this Hamiltonian has the exciton energies for the absorption for the individual dye molecule on the diagonal. The off-diagonal elements determine the coupling between the dye molecules, which in the dipole-dipole coupling approximation is proportional to $R^{-3}_{ij}$, where $R_{ij}$ is the distance between molecules i and j. The absorption spectrum is defined by the eigenvalues ($\epsilon_i$) and eigenvectors ($c_{ij}$) of this Hamiltonian:
I(E)=\sum_{i}\delta(\epsilon_i-E)\left|\sum_{j}c_{ij}\right|^2  
The delta function can be implemented by binning the intensity of the transitions to each eigenstate in a histogram depending on its eigen energy.
At GitHub.com a serial program doing this calculation can be found. 
1. Identify the subprocesses that can be made parallel. 
2. Implement OpenMP parallelisation
3. Test the scaling by running the program with different number of threads (If the calculation is very fast or very slow on your computer you can change the number of dye molecules in the aggregate to get a reasonable time duration of the calculations.)
4. Do you expect the scaling to be better or worse if you double the number of dye molecules? 
