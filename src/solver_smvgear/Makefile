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
	$(AR) crs $(LIBDIR)libsmvgear.a $(OBJ)

$(OBJDIR)%.o:%.f
	$(FC) $(FCFLAGS) $(MODFLAG) -c $< -o $@

clean:
	@rm -f $(MOD)
	@rm -f $(OBJ)
	@rm -f $(LIBDIR)libsmvgear.a

$(OBJDIR)mod_smvgear_comode.o : mod_smvgear_comode.f
$(OBJDIR)mod_smvgear_core.o : mod_smvgear_core.f mod_smvgear_comode.o
$(OBJDIR)mod_smvgear_interface.o : mod_smvgear_interface.f mod_smvgear_comode.o mod_smvgear_core.o	
