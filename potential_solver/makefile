FC=gfortran
FFLAGS=-O3 -Wall -Wextra -fdefault-real-8 -fcheck=all -fbackslash

SRC=constants.f90 inout.f90 pde_routines.f90 solver.f90
OBJ=${SRC:.f90=.o}

SRC2=constants.f90 inout.f90 pde_routines.f90 testing.f90
OBJ2=${SRC2:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

fdefault: solution test

solution: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

test: $(OBJ2)
	$(FC) $(FFLAGS) -o $@ $(OBJ2)

clean: $(OBJ)
	rm *.o *.mod *.dat

clean_all: $(OBJ)
	rm *.o *.mod *.dat solution
