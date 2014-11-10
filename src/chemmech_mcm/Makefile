OBJDIR := ../../obj/
MODDIR := ../../mod/
MODFLAG:= -J$(MODDIR)
LIBDIR := ../../lib/
include ../Makefile_Header

SRC := $(wildcard *.f)
SRC2 := $(wildcard *.f90)
OBJ := $(patsubst %.f,$(OBJDIR)%.o,$(SRC)) $(patsubst %.f90,$(OBJDIR)%.o,$(SRC2))
MOD := $(patsubst %.f,$(MODDIR)%.mod,$(SRC)) $(patsubst %.f90,$(MODDIR)%.mod,$(SRC2))

vpath %.o $(OBJDIR)

all:$(OBJ)
	$(AR) crs $(LIBDIR)libmcm.a $(OBJ)
#	$(FC) $(FCFLAGS) -o $(EXE) $(OBJ)
#	@mv $(EXE) $(EXEDIR)
#	@ln -s $(EXEDIR)$(EXE) .

$(OBJDIR)%.o:%.f
	$(FC) $(FCFLAGS) $(MODFLAG) -c $< -o $@

$(OBJDIR)%.o:%.f90
	$(FC) $(FC90FLAGS) $(MODFLAG) -c $< -o $@

clean:
	@rm -f $(MOD)
	@rm -f $(OBJ)
	@rm -f $(LIBDIR)libmcm.a

$(OBJDIR)module_mcm_cheminfo.o:module_mcm_cheminfo.f 
$(OBJDIR)module_mcm_rate.o : module_mcm_rate.f90 
