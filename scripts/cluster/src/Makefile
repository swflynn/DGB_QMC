#==============================================================================#
#			   Fortran Compiler/Flags
#-O0=no optimization,-O3=full optimization
#==============================================================================#
FC=gfortran 
FFLAGS=-O2
#FFLAGS=-O0 -Wall -Wextra -pedantic -fcheck=all -fbacktrace -llapack
#==============================================================================#
#			   	Targets 
#==============================================================================#
all: dgb_cluster
#==============================================================================#
#			Compile all Fortran Files
#==============================================================================#
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $^
#==============================================================================#
#				DGB Cluster
#requires mbpol directory/library to run their surface (mbpol/libmbpol.a)
#==============================================================================#
dgb_cluster: sobol.o sobol_stdnormal.o TIP4P.o mbpol/libmbpol.a dgb_cluster.o 
	$(FC) -o $@ $^ mbpol/libmbpol.a -lstdc++ -llapack
#==============================================================================#
clean:
	rm -rf *.o *.mod dgb_cluster
