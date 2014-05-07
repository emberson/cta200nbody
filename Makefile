FC=ifort
FFLAGS=-coarray

cube.x: cube.f90
	$(FC) $(FFLAGS) -coarray-num-images=8 $< -o $@

clean:
	rm -f *.o *.x

