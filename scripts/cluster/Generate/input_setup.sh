#=============================================================================80
#                           Generate Input Files 
#==============================================================================!
#    Discussion:
#Make a new directory (suggested name=Number Basis Functions)
#Generate various subdirectories with varying alpha parameter
#==============================================================================!
#   Modified:
# 4 April 2019
#   Author:
# Shane Flynn
#==============================================================================!
#NG         ==> Number of Gaussian Basis Functions
#potential  ==> Potential Surface
#Nsobol     ==> Number of sobol points for numerical integration
#freq       ==> Frequency at which the generalized eigenvalue problem is solved
#Ndir       ==> Total number of directories to generate
#i_alpha    ==> Initial Alpha0 value to start increment from
#inc        ==> increment to increase Alplha0 by 
#==============================================================================!
mkdir 10
cd 10
NG=10;
potential=tip4p;
Nsobol=200;
freq=10;
counter=1;
Ndir=5;
i_alpha=0.0;
inc=0.5;
while [ "$counter" -le "$Ndir" ]; do
  mkdir "$counter"
  cp ../cage.xyz "$counter"
  cp ../submit.sh "$counter"
  cd "$counter"
#==============================================================================!
#   single arrow makes a file, double arrow appends to an existing file
#==============================================================================!
  echo "$potential" > input
#==============================================================================!
#           can't use floating point with bash, use bc instead
#==============================================================================!
  echo "$i_alpha+ $counter*$inc" | bc  >> input
  echo "$NG" >> input
  echo "$Nsobol" >> input
  echo "$freq" >> input
  echo cage.xyz >> input
  echo 1 >> input
  cd ..
  counter=`expr "$counter" + 1`;
done
