#=============================================================================80
#                       Plot Eigenvalues vs Alpha0 
#==============================================================================!
#    Discussion:
#Xmgrace plot all n-th eigenvalues wrt alpha0 (1 Plot)
#==============================================================================!
#   Modified:
# 4 April 2019
#   Author:
# Shane Flynn
#==============================================================================!
#for i in {initial directory number..final directory number};
#-nxy plots all data assuming:: x y1  y2  ... y_n
#==============================================================================!
for i in {1..150}; do echo -n " -nxy $i/alpha.dat "; done | xargs xmgrace
