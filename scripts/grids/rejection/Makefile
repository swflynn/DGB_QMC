OBJS = Sobol S_std rej main
FC = gfortran -g -fcheck=all

all: Sobol S_std rej main

Sobol:
	$(FC) -c sobol.f90

S_std:
	$(FC) -c sobol_stdnormal.f90

rej:
	$(FC) -c rejection.f90 

main:
	$(FC) -O rejection.o sobol.o sobol_stdnormal.o 

clean:
	rm -rf *.o *.mod *.out $(OBJS)
