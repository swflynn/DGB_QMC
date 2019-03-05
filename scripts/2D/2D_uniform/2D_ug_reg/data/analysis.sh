# run the fortran code
./a.out < input.dat
# generate data analysis
python analysis.py eigenvalues.dat true.dat
# move everything to new directory for plotting
mv all.dat centers.dat eigenvalues.dat overlap_eigenvalues.dat overlap_regularized.dat simulation.dat true.dat reg
cd reg
head -21 all.dat > plot.dat
