OBJS = Sobol S_std 2D main
FC = gfortran -g -fcheck=all

all: Sobol S_std 2D main

Sobol:
	$(FC) -c sobol.f90

S_std:
	$(FC) -c sobol_stdnormal.f90

2D:
	$(FC) -c 2D_ug.f90 -llapack

main:
	$(FC) -O 2D_ug.o sobol.o sobol_stdnormal.o -llapack

clean:
	rm -rf *.o *.mod *.log *.out $(OBJS)
