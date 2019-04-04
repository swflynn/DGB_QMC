#=============================================================================80
#                           Generate Input Files 
#==============================================================================!
#    Discussion:
#Search all subdirectories, if they contain the cluster submission (submit.sh) 
#file then move to that directory and submit (qsub)
#==============================================================================!
#   Modified:
# 4 April 2019
#   Author:
# Shane Flynn
#==============================================================================!
for d in `find . -type d`
do ( cd "$d"
     if test ! -f submit.sh; then continue; fi
     qsub submit.sh 
) done
