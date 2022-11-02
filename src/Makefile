# These are good settings for debugging. Will produce lots of warnings etc
# but its very slow
#FC = ifort
#FF = -g -O0 -traceback -qopenmp -check bounds -checkformat -check stack -check uninit  -warn

# Try gfortran
FC = gfortran
FF = -g -fbacktrace -fopenmp -fcheck=all -O0 -ffpe-trap=invalid,zero,overflow
CPP = #-fpp

# TODO: Switch to more optimized code, e.g. O3 and remove checks and -g 

# If you have nf-config, this works
NC_LIB=$(shell nf-config --flibs)
NC_INC=$(shell nf-config --cflags)

# If not, set netcdf dirs yourself
#NC_LIB=-L /usr/lib/ -lnetcdf -lnetcdff
#NC_INC=-I /usr/include/

strf = -Dpsiyz -Dpsixz -Dpsixr -Dpsiyz -Dpsiyr -Dpsirr

objects := support.o vars.o grid.o data.o psi_main.o

%.o : %.F90 
	$(FC) $(FF) -c $(CPP) $(strf) $(NC_INC) $(NC_LIB) $< -o $@

hydrocode.x : $(objects)
	$(FC) $(FF) -o hydrocode.x $(objects) $(NC_INC) $(NC_LIB)

clean:
	rm -rf *.o *.mod
	rm -rf hydrocode.x 

