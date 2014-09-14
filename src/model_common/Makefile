OBJDIR := ../../obj/
MODDIR := ../../mod/

include ../Makefile_Header
	
SRC := $(wildcard *.f)
OBJ := $(patsubst %.f,$(OBJDIR)%.o,$(SRC))
MOD := $(patsubst %.f,$(MODDIR)%.mod,$(SRC))

vpath %.o $(OBJDIR)

all:$(OBJ)

$(OBJDIR)%.o:%.f
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	@rm $(MOD)
	@rm $(OBJ)

$(OBJDIR)module_chemmech_common.o : module_model_cheminfo.f module_model_parameter.o
$(OBJDIR)module_model_parameter.o : module_model_parameter.f
