OBJDIR := ../../obj/
MODDIR := ../../mod/
MODFLAG:= -J$(MODDIR)

include ../Makefile_Header
	
SRC := $(wildcard *.f)
OBJ := $(patsubst %.f,$(OBJDIR)%.o,$(SRC))
MOD := $(patsubst %.f,$(MODDIR)%.mod,$(SRC))

vpath %.o $(OBJDIR)

all:$(OBJ)

$(OBJDIR)%.o:%.f
	$(FC) $(FCFLAGS) $(MODFLAG) -c $< -o $@

clean:
	@rm -f $(MOD)
	@rm -f $(OBJ)

$(OBJDIR)module_chemmech_common.o : module_chemmech_common.f module_model_parameter.o
$(OBJDIR)module_conc_common.o : module_conc_common.f module_model_parameter.o
$(OBJDIR)module_domain_common.o:module_domain_common.f module_model_parameter.o
$(OBJDIR)module_time_common.o:module_time_common.f module_model_parameter.o julday_mod.o util_time_mod.o
$(OBJDIR)module_model_option.o :module_model_option.f
$(OBJDIR)module_model_parameter.o : module_model_parameter.f
