OBJDIR := ../../obj/
MODDIR := ../../mod/

include ../Makefile_Header

SRC := $(wildcard *.f)
OBJ := $(patsubst %.f,$(OBJDIR)%.o,$(SRC))
MOD := $(patsubst %.f,$(MODDIR)%.mod,$(SRC))

vpath %.o $(OBJDIR)

all:$(OBJ)
	@echo $(OBJ)

$(OBJDIR)%.o:%.f
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	@rm $(MOD)
	@rm $(OBJ)

$(OBJDIR)module_geoschem_rxntype.o : module_geoschem_rxntype.f module_model_parameter.o
$(OBJDIR)module_geoschem_cheminfo.o : module_geoschem_cheminfo.f module_geoschem_rxntype.o module_model_parameter.o
$(OBJDIR)module_geoschem_io.o : module_geoschem_io.f module_geoschem_cheminfo.o module_geoschem_rxntype.o module_model_parameter.o