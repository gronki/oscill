
FFLAGS=-mieee-fp -Wall -g

build: bin/oscill-base bin/oscill-one

bin/oscill-base: cgs_constants.o qho.o oscill-base.f90
		mkdir -p bin
		gfortran $(FFLAGS) -o $@ $^

bin/oscill-one: cgs_constants.o qho.o oscill-one.f90
		mkdir -p bin
		gfortran $(FFLAGS) -o $@ $^

cgs_constants.o: cgs_constants.f90
		gfortran $(FFLAGS) -o $@ -c $<

qho.o: qho.f90
		gfortran $(FFLAGS) -o $@ -c $<

clean:
		rm -f *.o *.mod
