# run the fortran code
./a.out < input
# generate data analysis
python analysis.py eigenvalues.dat true.dat
# move everything to new directory for plotting
mv all.dat centers.dat eigenvalues.dat overlap_eigenvalues.dat simulation.dat true.dat data
cd data
head -21 all.dat > plot.dat
