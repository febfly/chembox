OBJDIR := ../../obj/
MODDIR := ../../mod/
MODFLAG:= -J$(MODDIR)
BINDIR := ../../bin/
LIBDIR := ../../lib/
LIBFLAG:= -L$(LIBDIR) -lcommon -lutil -lgeoschem -lream -lmcm -lfastjx -lsmvgear
EXE    := chembox.exe

include ../Makefile_Header
	
SRC := $(wildcard *.f)
OBJ := $(patsubst %.f,$(OBJDIR)%.o,$(SRC))
MOD := $(patsubst %.f,$(MODDIR)%.mod,$(SRC))

vpath %.o $(OBJDIR)

all:$(OBJ)
	$(LD) $(OBJ) $(LIBFLAG) -o $(EXE)
	@cp $(EXE) $(BINDIR)

$(OBJDIR)%.o: %.f
	$(FC) $(FCFLAGS) $(MODFLAG) -c $< -o $@

clean:
	@rm -f $(MOD)
	@rm -f $(OBJ)
	@rm -f $(EXE)
	@rm -f $(BINDIR)$(EXE)

$(OBJDIR)chembox.o: chembox.f module_init.o module_chem.o
$(OBJDIR)module_init.o:module_init.f
$(OBJDIR)module_chem.o:module_chem.f
