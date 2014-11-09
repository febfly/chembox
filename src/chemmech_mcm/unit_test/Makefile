
EXE    := test_mcm.exe

include ../../Makefile_Header

SRC := $(wildcard *.f)
OBJ := $(patsubst %.f,%.o,$(SRC))
SRC2:= $(wildcard *.f90)
OBJ := $(OBJ) $(patsubst %.f90,%.o,$(SRC2))

#SRCUNIT := $(wildcard $(UNITDIR)*.f)
#SRCUNIT := $(filter-out $(UNITDIR)module_ream_common.f,$(SRCUNIT))
#OBJUNIT := $(notdir $(patsubst %.f,%.o,$(SRCUNIT)))


all:$(OBJ)
	@echo unit test!
	$(LD)  $(OBJ) $(LINK) -o $(EXE)

%.o:%.f
	$(FC) $(FCFLAGS) -c $< -o $@

%.o:%.f90
	$(FC) $(FC90FLAGS) -c $< -o $@

clean:
	@rm *.mod
	@rm *.o

drive.o:drive.f module_mcm_cheminfo.o module_mcm_rate.o
module_mcm_cheminfo.o:module_mcm_cheminfo.f
module_mcm_rate.o:module_mcm_rate.f90

