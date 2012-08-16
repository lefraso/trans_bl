
OBJ1 = constants.f90 solvers.f90 derivative.f90 at_wall.f90 coef_generation.f90 

OBJ2 = constants.f90 beta_gen.f90

coef : $(OBJ1)
	gfortran $(OBJ1) -o coef

beta : $(OBJ2)
	gfortran $(OBJ2) -o beta

clean :
	rm *.bin *.mod *~ coef beta fs *.dat

fs:	
	./fs.sh

.f.o:
	gfortran -c -02 $<
