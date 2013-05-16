.PHONY : all clean

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

EXE=GearProc.exe

F90FLAGS := $(F90FLAGS) -DDOREORD=$(DOREORD) -DBLOCKSIZE=$(BLOCKSIZE)
all: $(OBJS) $(EXE)

%.o : %.F90
	$(FC) $(F90FLAGS) -c $<

$(EXE): $(OBJS)
	$(FC) $(F90FLAGS) $(PFUNIT_LFLAGS) -o $@ $(OBJS)
clean:
	rm -f *.o *.mod *.x *.exe
