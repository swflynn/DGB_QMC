# write the first n eigenvalues to a new file for plotting
# for water only care about bend and stretch ~ 3000cm-1 from ground state
for d in `find . -type d`
do ( cd "$d"
     if test ! -f alpha_conv.dat; then continue; fi
    awk -v n=100 'n==c{exit}n-c>=NF{print;c+=NF;next}{for(i=1;i<=n-c;i++)printf "%s ",$i;print x;exit}' alpha_conv.dat > alpha.dat
) done
