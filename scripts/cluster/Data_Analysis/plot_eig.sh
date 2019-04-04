#=============================================================================80
#                       Plot Eigenvalues vs Nsobol
#==============================================================================!
#    Discussion:
#Xmgrace plot all n-th eigenvalues wrt Nsobol Iteration (for a single Alpha0)
#==============================================================================!
#   Modified:
# 4 April 2019
#   Author:
# Shane Flynn
#==============================================================================!
#eigs.dat   ==> Nsobol-iteration E0 E1  ... En
#-nxy plots all data assuming:: x y1  y2  ... yn columns
#==============================================================================!
for d in `find . -type d`
do ( cd "$d"
     if test ! -f eigs.dat; then continue; fi
     xmgrace -batch ../eig.bat -nxy eigs.dat -nosafe -hardcopy
) done
