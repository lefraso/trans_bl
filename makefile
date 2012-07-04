# makefile for gvts
OBJ1 = main_par.o derivs_par.o loop_par.o nlterms_par.o escreve_par.o penta_par.o poisson_par.o filter_par.o immersed_par.o
OBJ2 = integ.o
OBJ3 = ftanalisys_par.o
OBJ4 = isoq_par.o

ru : $(OBJ1)
	mpif77 $(OBJ1) -o rugos

int : $(OBJ2)
	mpif77 $(OBJ2) -o integ

ft : $(OBJ3)
	mpif77 $(OBJ3) -o ft

iso : $(OBJ4)
	mpif77 $(OBJ4) -o isoq

.f.o:
	mpif77 -c -O3 $<
#	mpif77 -O0 -g -Wall -fbounds-check -c $<

clean:
	rm *~ *.o *bin *dat fs rugos
