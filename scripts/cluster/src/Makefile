OBJS = Tip Sobol S_std dgb main
FC = gfortran 

all: Tip Sobol S_std dgb main

Tip:
	$(FC) -c  TIP4P.f90

Sobol:
	$(FC) -c  sobol.f90

S_std:
	$(FC) -c  sobol_stdnormal.f90

dgb:
	$(FC) -c -O dgb.f90  -llapack
#	$(FC) -c -g -O0 -Wall -pedantic -fcheck=all -fbacktrace dgb.f90  -llapack

main:
	$(FC) -O dgb.o  TIP4P.o sobol.o sobol_stdnormal.o mbpol/libmbpol.a -lstdc++ -llapack
#	$(FC) -O -Wall -fcheck=all -g -fbacktrace dgb.o TIP4P.o sobol.o sobol_stdnormal.o mbpol/libmbpol.a -lstdc++ -llapack

clean:
	rm -rf *.o *.mod *.out $(OBJS)
