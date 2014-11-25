# target
TAR  = square
OBJ = square.o

# memory increase on mac
ifeq ($(OSX),true)
STACK = -Wl,-stack_size -Wl,0x40000000
else
STACK = 
endif 

# double precision flag
ifeq ($(DUB),true)
PREC = -DDUBON -fdefault-real-8
else
PREC = 
endif 

# debugging flags
#FLAGS = -g -Wall 
# production code flags
FLAGS =  $(PREC) -Ofast -funroll-loops --param max-unroll-times=4 $(STACK)

# MPI flag
ifeq ($(MPI),true)
FC = mpif90 -cpp -DMPION
else
FC = $(F90) -cpp
endif

# library info
LIBS = -lfluid -lgeom
LDIR = -L$(MGLHOME)/oop/. -L$(MGLHOME)/geom_lib/.
IDIR = -I$(MGLHOME)/oop/fluid_lib_include -I$(MGLHOME)/geom_lib/geom_lib_include

# make
.SUFFIXES: .f90
.f90.o:
	$(FC) $(FLAGS) $(IDIR) -c $*.f90

$(TAR): $(OBJ)
	$(FC) $(FLAGS) -o $(TAR) $(OBJ) $(LDIR) $(LIBS)

all: $(TAR)

clean: 
	rm -f *.o *~ *.mod

cleaner:
	rm -f *.o *~ *.mod $(TAR) *.nc fort.*

