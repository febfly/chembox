OBJDIR := ../../obj/
MODDIR := ../../mod/
MODFLAG:= -J$(MODDIR)
LIBDIR := ../../lib/
include ../Makefile_Header

SRC := $(wildcard *.f)
OBJ := $(patsubst %.f,$(OBJDIR)%.o,$(SRC))
MOD := $(patsubst %.f,$(MODDIR)%.mod,$(SRC))

vpath %.o $(OBJDIR)

all:$(OBJ)
	$(AR) crs $(LIBDIR)libream.a $(OBJ)
#	$(FC) $(FCFLAGS) -o $(EXE) $(OBJ)
#	@mv $(EXE) $(EXEDIR)
#	@ln -s $(EXEDIR)$(EXE) .

$(OBJDIR)%.o:%.f
	$(FC) $(FCFLAGS) $(MODFLAG) -c $< -o $@

clean:
	@rm -f $(MOD)
	@rm -f $(OBJ)
	@rm -f $(LIBDIR)libream.a

$(OBJDIR)module_ream_common.o:module_ream_common.f module_chemmech_common.o
$(OBJDIR)module_ream_rxntype.o : module_ream_rxntype.f module_ream_common.o  module_chemmech_common.o 
$(OBJDIR)module_ream_cheminfo.o : module_ream_cheminfo.f module_ream_rxntype.o module_ream_common.o  module_chemmech_common.o
$(OBJDIR)module_ream_io.o : module_ream_io.f module_ream_cheminfo.o module_ream_rxntype.o module_ream_common.o  module_chemmech_common.o
