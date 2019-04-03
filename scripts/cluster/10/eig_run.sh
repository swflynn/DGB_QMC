for d in `find . -type d`
do ( cd "$d"
     if test ! -f eigs.dat; then continue; fi
     xmgrace -batch ../eig.bat -nxy eigs.dat -nosafe -hardcopy
) done
