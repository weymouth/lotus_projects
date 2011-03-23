#
# makefile for BDIM
#
ifeq ($(MPI),true)
VAL = par
else
VAL = ser
endif

MGL = $(MGLHOME)

TAR1 = exec_BDIM
OBJ1 = fglo.o f$(VAL).o fgri.o fino.o fana.o figs.o fnur.o fgeo.o fran.o fslp.o fbod.o fblo.o fvof.o flin.o fqui.o fdom.o fmul.o fpre.o fapr.o main.o

TAR2 = exec_grid
OBJ2 = design_grid.o

TAR3 = exec_geom
OBJ3 = fana.o design_geom.o

TAR4 = exec_IGS
OBJ4 = fana.o figs.o fnur.o design_IGS.o

TAR5 = exec_apriori
OBJ5 = fblo.o apriori.o

FLAGS = -O3 -Wl,-stack_size,0x10000000

ifeq ($(F90),gfortran)
F90FLAGS = -fdefault-real-8 $(FLAGS)
endif

ifeq ($(F90),ftn)
F90FLAGS = -r8 -O3
else
F90 = ifort
F90FLAGS = -r8 $(FLAGS)
endif

ifeq ($(MPI),true)
FC = mpif90 -f90=$(F90)
else
FC = $(F90)
endif

ifeq ($(F90),ftn)
FC = ftn
endif

# library info
LIBS = -lmitmg_$(VAL)
LDIR = -L$(MGLHOME)/mglibrary/.
IDIR = -I$(MGLHOME)/mglibrary/mitmg_$(VAL)_include

#
.SUFFIXES: .f90
.f90.o:
	$(FC) $(F90FLAGS) $(IDIR) -c $*.f90

$(TAR1): $(OBJ1)
	$(FC) $(F90FLAGS) -o $(TAR1) $(OBJ1) $(LDIR) $(LIBS)

pre: $(TAR2) $(TAR3) $(TAR4)

post: $(TAR5)

$(TAR2): $(OBJ2)
	$(F90) $(F90FLAGS) -o $(TAR2) $(OBJ2)

$(TAR3): $(OBJ3)
	$(F90) $(F90FLAGS) -o $(TAR3) $(OBJ3)

$(TAR4): $(OBJ4)
	$(F90) $(F90FLAGS) -o $(TAR4) $(OBJ4)

$(TAR5): $(OBJ5)
	$(F90) $(F90FLAGS) -o $(TAR5) $(OBJ5)


clean:
	-rm -f *.o 
	-rm -f *~
	-rm -f *.mod

fullclean:
	-rm -f *.o 
	-rm -f *~
	-rm -f *.mod
	-rm -f exec_*
