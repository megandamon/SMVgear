.PHONY : all gearproc clean

SRCS=$(wildcard *.F90)
OBJS=$(SRCS:.F90=.o)
FC=ifort

BLOCKSIZE=32
DOOPT = 1
DOREORD = 0
PROFILE = 0

F90FLAGS := -O3 -vec-report0 -ftz -align all -fno-alias -convert big_endian -fPIC -fpe0 -fp-model precise  -align dcommons -check bounds -check uninit
#F90FLAGS := -O0 -g -traceback -vec-report0 -ftz -align all -fno-alias -convert big_endian -fPIC -fpe0 -fp-model precise  -align dcommons
#F90FLAGS := -O3 -vec-report0 -ftz -align all -fno-alias -convert big_endian -fPIC -fpe0 -fp-model precise  -align dcommons -check uninit

PFUNIT=/Users/mrdamon/git/pfunit_install/
PFUNIT_LFLAGS += -L$(PFUNIT)/lib -lpfunit


EXE=Do_Smv2_Solver.exe

F90FLAGS := $(F90FLAGS) -DDOREORD=$(DOREORD) -DBLOCKSIZE=$(BLOCKSIZE)
all: $(OBJS) $(EXE) gearproc

%.o : %.F90
	$(FC) $(F90FLAGS) -c $<

$(EXE): $(OBJS)
	$(FC) $(F90FLAGS) $(PFUNIT_LFLAGS) -o $@ $(OBJS)

gearproc:
	make -C GearProc all
	
smv2chem_solver.o:timing_mod.o SmvChem_mod.o
physproc.o:timing_mod.o GmiPrintError_mod.o
smvgear.o:GmiPrintError_mod.o GmiMechanism_mod.o GmiSparseMatrix_mod.o
GmiManager_mod.o:GmiPrintError_mod.o Smv2Chem2_mod.o
SmvChem_mod.o:GmiPrintError_mod.o
ChemTable_mod.o:Smv2Chem2_mod.o Smv2Chem1_mod.o

#timing_mod.o:timing_mod.F90
#	$(FC) -c $(F90FLAGS_R8) timing_mod.F90

clean:
	rm -f *.o *.mod *.x *.exe
	make -C GearProc clean

ifdef PFUNIT
  include $(PFUNIT)/include/pFUnit.makefile
endif
