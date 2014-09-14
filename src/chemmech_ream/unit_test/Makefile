UNITDIR := ../
EXE    := test_ream.exe

include ../../Makefile_Header

SRC := $(wildcard *.f)
OBJ := $(patsubst %.f,%.o,$(SRC))

SRCUNIT := $(wildcard $(UNITDIR)*.f)
SRCUNIT := $(filter-out $(UNITDIR)module_ream_common.f,$(SRCUNIT))
OBJUNIT := $(notdir $(patsubst %.f,%.o,$(SRCUNIT)))

vpath %.f $(UNITDIR)

all:$(OBJUNIT) $(OBJ)
	@echo unit test!
	$(LD) $(OBJUNIT) $(OBJ) $(LINK) -o $(EXE)

%.o:%.f
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	@rm *.mod
	@rm *.o

main.o:main.f module_ream_common.o module_ream_rxntype.o module_ream_cheminfo.o module_chemmech_common.o
module_chemmech_common.o:module_chemmech_common.f
module_ream_common.o:module_ream_common.f module_chemmech_common.o
module_ream_rxntype.o : module_ream_rxntype.f module_ream_common.o module_chemmech_common.o
module_ream_cheminfo.o : module_ream_cheminfo.f module_ream_rxntype.o module_ream_common.o module_chemmech_common.o
module_ream_io.o : module_ream_io.f module_ream_cheminfo.o module_ream_rxntype.o module_ream_common.o module_chemmech_common.o
