# run the fortran code
./a.out < input
# generate data analysis
python iter_err.py eigenvalues.dat theory.dat
mkdir data
# move everything to new directory for plotting
mv *.dat data
