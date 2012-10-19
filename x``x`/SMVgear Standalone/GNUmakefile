SRCS=$(wildcard *.F90)
OBJS=$(SRCS:.F90=.o)
FC=ifort

F90FLAGS := -O0 -traceback -vec-report0 -ftz -align all -fno-alias -convert big_endian -fPIC -fpe0 -fp-model precise  -align dcommons
F90FLAGS := -O3 -vec-report0 -ftz -align all -fno-alias -convert big_endian -fPIC -fpe0 -fp-model precise  -align dcommons

EXE=Do_Smv2_Solver.exe

all: $(OBJS) $(EXE)

%.o : %.F90
	$(FC) $(F90FLAGS) -c $<

$(EXE): $(OBJS)
	$(FC) $(F90FLAGS) -o $@ $(OBJS)

smv2chem_solver.o:timing_mod.o
physproc.o:timing_mod.o GmiPrintError_mod.o
smvgear.o:GmiPrintError_mod.o

#timing_mod.o:timing_mod.F90
#	$(FC) -c $(F90FLAGS_R8) timing_mod.F90

clean:
	rm -f *.o *.mod *.x
