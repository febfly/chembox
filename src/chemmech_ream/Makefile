OBJDIR := ../../obj/
MODDIR := ../../mod/

include ../Makefile_Header

SRC := $(wildcard *.f)
OBJ := $(patsubst %.f,$(OBJDIR)%.o,$(SRC))
MOD := $(patsubst %.f,$(MODDIR)%.mod,$(SRC))

vpath %.o $(OBJDIR)

all:$(OBJ)
	@echo $(OBJ)
#	$(FC) $(FCFLAGS) -o $(EXE) $(OBJ)
#	@mv $(EXE) $(EXEDIR)
#	@ln -s $(EXEDIR)$(EXE) .

$(OBJDIR)%.o:%.f
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	@rm $(MOD)
	@rm $(OBJ)

$(OBJDIR)module_model_parameter.o : module_model_parameter.f
$(OBJDIR)module_ream_rxntype.o : module_ream_rxntype.f module_model_parameter.o
$(OBJDIR)module_ream_cheminfo.o : module_ream_cheminfo.f module_ream_rxntype.o module_model_parameter.o
$(OBJDIR)module_ream_io.o : module_ream_io.f module_ream_cheminfo.o module_ream_rxntype.o module_model_parameter.o
