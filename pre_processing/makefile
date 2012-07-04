
OBJ1 = constants.f90 solvers.f90 derivative.f90 at_wall.f90 coef_generation.f90 

coef : $(OBJ1)
	gfortran $(OBJ1) -o coef

clean :
	rm *bin *mod *~ coef fs *dat

fs:	
	./fs.sh

.f.o:
	gfortran -c -02 $<
