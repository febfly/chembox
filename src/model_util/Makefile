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
	$(AR) crs $(LIBDIR)libutil.a $(OBJ)

$(OBJDIR)%.o:%.f
	$(FC) $(FCFLAGS) $(MODFLAG) -c $< -o $@

clean:
	@rm -f $(MOD)
	@rm -f $(OBJ)
	@rm -f $(LIBDIR)libutil.a

$(OBJDIR)julday_mod.o : julday_mod.f
$(OBJDIR)util_time_mod.o:util_time_mod.f julday_mod.o
$(OBJDIR)module_util.o: module_util.f
