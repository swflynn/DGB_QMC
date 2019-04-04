#=============================================================================80
#                       Convergence wrt. Sobol Points
#==============================================================================!
#    Discussion:
#Generally only care about the t lowest eigenvalues (not all NG eigenvalues) 
#For water ~3000cm-1 above the ground state corresponds to vibrations
#uses xmgrace to plot n eigenvalues for each alpha0 as a function of N-Sobol
#==============================================================================!
#   Modified:
# 4 April 2019
#   Author:
# Shane Flynn
#==============================================================================!
#eigenvalues.dat    ==> Nsobol  E0  E1  E2  .... Et  ......   NG
#eigs.dat           ==> Nsobol  E0  E1  E2  .... Et
#t                  ==> remove all eigenvalues from (t+1) to NG
#==============================================================================!
for d in `find . -type d`
do ( cd "$d"
    if test ! -f eigenvalues.dat; then continue; fi
    awk -v f=1 -v t=100 '{for(i=f;i<=t;i++) printf("%s%s",$i,(i==t)?"\n":OFS)}' eigenvalues.dat > eigs.dat
    ) done
