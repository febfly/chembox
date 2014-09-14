UNITDIR := ../
EXE    := test_geoschem.exe

include ../../Makefile_Header

SRC := $(wildcard *.f)
OBJ := $(patsubst %.f,%.o,$(SRC))

SRCUNIT := $(wildcard $(UNITDIR)*.f)
SRCUNIT := $(filter-out $(UNITDIR)module_geoschem_common.f,$(SRCUNIT))
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

main.o:main.f module_geoschem_common.o module_geoschem_rxntype.o module_geoschem_cheminfo.o module_chemmech_common.o
module_chemmech_common.o:module_chemmech_common.f
module_geoschem_common.o:module_geoschem_common.f module_chemmech_common.o
module_geoschem_rxntype.o : module_geoschem_rxntype.f module_geoschem_common.o module_chemmech_common.o
module_geoschem_cheminfo.o : module_geoschem_cheminfo.f module_geoschem_rxntype.o module_geoschem_common.o module_chemmech_common.o
module_geoschem_io.o : module_geoschem_io.f module_geoschem_cheminfo.o module_geoschem_rxntype.o module_geoschem_common.o module_chemmech_common.o
