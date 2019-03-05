# run the fortran code
./a.out < input.dat
# generate data analysis
#python iter_err.py eigenvalues.dat theory.dat
#mkdir data
# move everything to new directory for plotting
#mv error.dat centers.dat eigenvalues.dat overlap_eigenvalues.dat simulation.dat theory.dat data
#cd data
# generate file for plotting first 20 eigenvalues (for convenience)
#head -21 error.dat > plot_error.dat
