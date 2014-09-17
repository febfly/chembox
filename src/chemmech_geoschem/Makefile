OBJDIR := ../../obj/
MODDIR := ../../mod/
LIBDIR := ../../lib/
MODFLAG:= -J$(MODDIR)

include ../Makefile_Header

SRC := $(wildcard *.f)
OBJ := $(patsubst %.f,$(OBJDIR)%.o,$(SRC))
MOD := $(patsubst %.f,$(MODDIR)%.mod,$(SRC))

vpath %.o $(OBJDIR)

all:$(OBJ)
	$(AR) crs $(LIBDIR)libgeoschem.a $(OBJ)

$(OBJDIR)%.o:%.f
	$(FC) $(FCFLAGS) $(MODFLAG) -c $< -o $@

clean:
	@rm -f $(MOD)
	@rm -f $(OBJ)
	@rm -f $(LIBDIR)libgeoschem.a

$(OBJDIR)module_geoschem_common.o:module_geoschem_common.f module_chemmech_common.o
$(OBJDIR)module_geoschem_rxntype.o : module_geoschem_rxntype.f module_geoschem_common.o module_chemmech_common.o
$(OBJDIR)module_geoschem_cheminfo.o : module_geoschem_cheminfo.f module_geoschem_rxntype.o module_geoschem_common.o module_chemmech_common.o
$(OBJDIR)module_geoschem_io.o : module_geoschem_io.f module_geoschem_cheminfo.o module_geoschem_rxntype.o module_geoschem_common.o module_chemmech_common.o
