FC=gfortran
FFLAGS=-O3 -Wall -Wextra -fdefault-real-8 -fcheck=all

SRC=constants.f90 inout.f90 pde_routines.f90 solver.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

solution: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean: $(OBJ)
	rm *.o *.mod *.dat

clean_all: $(OBJ)
	rm *.o *.mod *.dat solution
