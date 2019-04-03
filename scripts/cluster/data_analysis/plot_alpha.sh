# plot all alpha.dat files
# xmgrace -param test.par -nxy 30/alpha.dat -nxy 40/alpha.dat 
# xmgrace -param test.par -nzy */alpha.dat 
for i in {1..150}; do echo -n " -nxy $i/alpha.dat "; done | xargs xmgrace
