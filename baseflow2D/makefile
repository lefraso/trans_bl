# makefile for baseflow
OBJ1 = main_par_2D.o derivs_par_2D.o loop_par_2D.o penta_par_2D.o poisson_par_2D.o
OBJ2 = integ.o

base : $(OBJ1)
	mpif77 $(OBJ1) -o basef

int : $(OBJ2)
	mpif77 $(OBJ2) -o integ

.f.o:
	mpif77 -c -O2 $<
#	mpif77 -O0 -g -Wall -fbounds-check -c $<

clean:
	rm *~ *.o *.bin *.dat basef
