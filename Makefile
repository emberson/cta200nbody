FC=ifort
FFLAGS=-mkl -coarray -coarray-num-images=8

cube.x: cube.f90 mkl_fftvec.o
	$(FC) $(FFLAGS) cube.f90 mkl_fftvec.o -o $@

mkl_fftvec.o: mkl_fftvec.f90
	$(FC) -fpp -c $<

clean:
	rm -f *.o *.x

