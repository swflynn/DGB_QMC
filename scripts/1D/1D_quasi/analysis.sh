python data.py eigenvalues.dat theory.dat
mv *.dat data
cp input data
cd data
xmgrace -param 1d_par.par centers.dat
xmgrace -settype xydy -param 1D_final.par final.dat
