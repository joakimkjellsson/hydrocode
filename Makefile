SYS    = vagn
SRC    = ./src/
NC_INC = $(shell nf-config --fflags)
NC_LIB = $(shell nf-config --flibs) 
PREP   = -Dpsiyr -Dpsiyz -Dpsirr 

ifeq ($(SYS),mac)
	F90  = gfortran
	OPT  = -O2
	OMP  = -fopenmp 
endif
ifeq ($(SYS),vagn)
	F90  = ifort
	OPT  = -O2 -funroll-loops -traceback -g
	OMP  = -openmp
	NC_LIB = -L/software/apps/netcdf/4.1.2/i1203/lib/ -lnetcdf -lnetcdff
	NC_INC = -I/software/apps/netcdf/4.1.2/i1203/include/
endif
ifeq ($(SYS),trio)
	F90  = ifort
	OPT  = -O2 -funroll-loops -xAVX
	OMP  = -openmp
endif

FLAGS = $(OPT) $(OMP) $(NC_INC) $(NC_LIB) $(PREP)

objects = vars.o grid.o support.o data.o psi_main.o

psi.x : $(objects)
	$(F90) -o psi.x $(LIB) $(INC) $(FLAGS) $(objects)

support.o : $(SRC)/support.f90
	$(F90) -fpp -c $(LIB) $(INC) $(FLAGS) $(SRC)/support.f90

vars.o : $(SRC)/vars.f90
	$(F90) -fpp -c $(LIB) $(INC) $(FLAGS) $(SRC)/vars.f90

data.o : $(SRC)/data.f90
	$(F90) -fpp -c $(LIB) $(INC) $(FLAGS) $(SRC)/data.f90

grid.o : $(SRC)/grid.f90
	$(F90) -fpp -c $(LIB) $(INC) $(FLAGS) $(SRC)/grid.f90

psi_main.o : $(SRC)/psi_main.f90
	$(F90) -fpp -c $(LIB) $(INC) $(FLAGS) $(SRC)/psi_main.f90

clean : 
	rm -rf *.o
	rm -rf *.mod
	rm psi.x
