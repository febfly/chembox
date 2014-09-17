OBJDIR := ../../obj/
MODDIR := ../../mod/
MODFLAG:= -J$(MODDIR)
LIBDIR := ../../lib/

include ../Makefile_Header

SRC := $(wildcard *.f90)
OBJ := $(patsubst %.f90,$(OBJDIR)%.o,$(SRC))
MOD := $(patsubst %.f90,$(MODDIR)%.mod,$(SRC))

vpath %.o $(OBJDIR)

all:$(OBJ)
	$(AR) crs $(LIBDIR)libfastjx.a $(OBJ)

$(OBJDIR)%.o:%.f90
	$(FC) $(FCFLAGS) $(MODFLAG) -c $< -o $@

clean:
	@rm -f $(MOD)
	@rm -f $(OBJ)
	@rm -f $(LIBDIR)libfastjx.a

$(OBJDIR)fjx_cmn_mod.o : fjx_cmn_mod.f90
$(OBJDIR)fjx_sub_mod.o : fjx_sub_mod.f90 fjx_cmn_mod.o
$(OBJDIR)fjx_init_mod.o : fjx_init_mod.f90 fjx_cmn_mod.o fjx_sub_mod.o

