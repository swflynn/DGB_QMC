#=============================================================================80
#                       Convergence wrt. Alpha0
#==============================================================================!
#    Discussion:
#Generally only care about the n lowest eigenvalues (not all NG eigenvalues) 
#For water ~3000cm-1 above the ground state corresponds to vibrations
#uses xmgrace to plot n eigenvalues for each alpha0
#==============================================================================!
#   Modified:
# 4 April 2019
#   Author:
# Shane Flynn
#==============================================================================!
#alpha_conv.dat     ==> alpha0 E0  E1  E2  .... En  ......   NG
#alpha.dat          ==> alpha0 E0  E1  E2  .... En
#==============================================================================!
for d in `find . -type d`
do ( cd "$d"
    if test ! -f alpha_conv.dat; then continue; fi
    awk -v n=100 'n==c{exit}n-c>=NF{print;c+=NF;next}{for(i=1;i<=n-c;i++)printf "%s ",$i;print x;exit}' alpha_conv.dat > alpha.dat
) done
