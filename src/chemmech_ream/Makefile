OBJDIR := ../../obj/
MODDIR := ../../mod/
MODFLAG:= -J$(MODDIR)
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
	$(FC) $(FCFLAGS) $(MODFLAG) -c $< -o $@

clean:
	@rm $(MOD)
	@rm $(OBJ)

$(OBJDIR)module_ream_parameter.o:module_ream_parameter.f
$(OBJDIR)module_ream_rxntype.o : module_ream_rxntype.f module_ream_parameter.o
$(OBJDIR)module_ream_cheminfo.o : module_ream_cheminfo.f module_ream_rxntype.o module_ream_parameter.o
$(OBJDIR)module_ream_io.o : module_ream_io.f module_ream_cheminfo.o module_ream_rxntype.o module_ream_parameter.o
