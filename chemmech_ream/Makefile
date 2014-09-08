include ../Makefile_Header
	
SRC := $(wildcard *.f)
OBJ := $(SRC:.f=.o)

all:$(OBJ)
#	$(FC) $(FCFLAGS) -o $(EXE) $(OBJ)
#	@mv $(EXE) $(EXEDIR)
#	@ln -s $(EXEDIR)$(EXE) .

clean:
	@rm *.o

module_model_parameter.o : module_model_parameter.f
module_ream_rxntype.o : module_ream_rxntype.f module_model_parameter.o
module_ream_cheminfo.o : module_ream_cheminfo.f module_ream_rxntype.o module_model_parameter.o
module_ream_io.o : module_ream_io.f module_ream_cheminfo.o module_ream_rxntype.o module_model_parameter.o
