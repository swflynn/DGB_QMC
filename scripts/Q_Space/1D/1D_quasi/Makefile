OBJS = Sobol S_std 1D main
FC = gfortran -g -fcheck=all

all: Sobol S_std 1D main

Sobol:
	$(FC) -c sobol.f90

S_std:
	$(FC) -c sobol_stdnormal.f90

1D:
	$(FC) -c 1D_qg.f90 -llapack

main:
	$(FC) -O 1D_qg.o sobol.o sobol_stdnormal.o -llapack

clean:
	rm -rf *.o *.mod *.out $(OBJS)
