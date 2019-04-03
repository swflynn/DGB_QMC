# write the first n eigenvalues to a new file for plotting
# for water only care about bend and stretch ~ 3000cm-1 from ground state
for d in `find . -type d`
do ( cd "$d"
    if test ! -f eigenvalues.dat; then continue; fi
    awk -v f=1 -v t=100 '{for(i=f;i<=t;i++) printf("%s%s",$i,(i==t)?"\n":OFS)}' eigenvalues.dat > eigs.dat
    ) done
