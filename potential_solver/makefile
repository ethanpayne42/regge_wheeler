FC=gfortran
FFLAGS=-O3 -Wall -Wextra -fdefault-real-8

SRC=constants.f90 inout.f90 pde_routines.f90 solver.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

solution: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean: $(OBJ)
	rm *.o *.mod

clean_all: $(OBJ)
	rm *.o *.mod solution
