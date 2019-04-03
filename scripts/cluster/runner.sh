for d in `find . -type d`
do ( cd "$d"
     if test ! -f submit.sh; then continue; fi
     qsub submit.sh 
) done
