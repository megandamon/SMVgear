SRCS=$(wildcard *.F90)
OBJS=$(SRCS:.F90=.o)
FC=ifort

F90FLAGS := -O0 -g -traceback -vec-report0 -ftz -align all -fno-alias -convert big_endian -fPIC -fpe0 -fp-model precise  -align dcommons
#F90FLAGS := -O3 -vec-report0 -ftz -align all -fno-alias -convert big_endian -fPIC -fpe0 -fp-model precise  -align dcommons -check bounds -check uninit
#F90FLAGS := -O3 -vec-report0 -ftz -align all -fno-alias -convert big_endian -fPIC -fpe0 -fp-model precise  -align dcommons -check uninit
PFUNIT_LFLAGS += -L$(PFUNIT)/lib -lpfunit


EXE=Do_Smv2_Solver.exe

all: $(OBJS) $(EXE)

%.o : %.F90
	$(FC) $(F90FLAGS) -c $<

$(EXE): $(OBJS)
	$(FC) $(F90FLAGS) $(PFUNIT_LFLAGS) -o $@ $(OBJS)

smv2chem_solver.o:timing_mod.o SmvChem_mod.o
physproc.o:timing_mod.o GmiPrintError_mod.o
smvgear.o:GmiPrintError_mod.o GmiMechanism_mod.o GmiSparseMatrix_mod.o
GmiManager_mod.o:GmiPrintError_mod.o
SmvChem_mod.o:GmiPrintError_mod.o

#timing_mod.o:timing_mod.F90
#	$(FC) -c $(F90FLAGS_R8) timing_mod.F90

clean:
	rm -f *.o *.mod *.x

ifdef PFUNIT
  include $(PFUNIT)/include/pFUnit.makefile
endif
