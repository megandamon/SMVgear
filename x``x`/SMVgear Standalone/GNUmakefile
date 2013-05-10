SRCS=$(wildcard *.F90)
OBJS=$(SRCS:.F90=.o)
FC=ifort

BLOCKSIZE=32
DOOPT = 1
DOREORD = 0
PROFILE = 0

# look at -O0, -g, etc.
#F90FLAGS := -O3 -vec-report0 -ftz -align all -fno-alias -convert big_endian -fPIC -fpe0 -fp-model precise  -align dcommons -check uninit

F90FLAGS := -openmp -O3 -fp-model source#-check all #-debug all
EXE=smv2solver.x

ifeq "$(DOOPT)" "1"
	F90FLAGS := $(F90FLAGS) -DOPT
endif
ifeq "$(PROFILE)" "1"
	F90FLAGS := $(F90FLAGS) -no-offload -profile-functions -profile-loops=all -profile-loops-report=2
endif

F90FLAGS := $(F90FLAGS) -DDOREORD=$(DOREORD) -DBLOCKSIZE=$(BLOCKSIZE) 

F90FLAGS := $(F90FLAGS) -align array64byte -align zcommons -align rec32byte
all: $(OBJS) $(EXE)

%.o : %.F90
	$(FC) $(F90FLAGS) -c $<

$(EXE): $(OBJS)
	$(FC) $(F90FLAGS) -o $@ $(OBJS)

smv2chem_solver.o:timing_mod.o 
physproc.o:timing_mod.o GmiPrintError_mod.o
smvgear.o:GmiPrintError_mod.o GmiMechanism_mod.o GmiSparseMatrix_mod.o
GmiManager_mod.o:GmiPrintError_mod.o

#timing_mod.o:timing_mod.F90
#	$(FC) -c $(F90FLAGS_R8) timing_mod.F90

clean:
	rm -f *.o *.mod *.x *_exit.proc0017
