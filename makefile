# makefile for gvts
OBJ1 = main_par.o derivs_par.o loop_par.o nlterms_par.o escreve_par.o penta_par.o poisson_par.o filter_par.o immersed_par.o les_par.o
OBJ2 = ftanalysis_par.o
OBJ3 = isoq_par.o
OBJ4 = formatted_par.o

code : $(OBJ1)
	mpif77 $(OBJ1) -o prog

ft : $(OBJ2)
	gfortran $(OBJ2) -o ft

iso : $(OBJ3)
	mpif77 $(OBJ3) -o isoq

form : $(OBJ4)
	gfortran $(OBJ4) -o form

.f.o:
	mpif77 -c -O3 -mcmodel=medium $<
#	mpif77 -O0 -g -Wall -fbounds-check -c $<

clean:
	rm *~ *.o *.bin *.dat fs prog ft isoq form saida pert_*
