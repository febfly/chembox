include ../Makefile_Header
OUTPATH :=
EXE := smvgear.exe
EXEDIR := 
	
SRC := $(wildcard *.f)
OBJ := $(SRC:.f=.o)

all:$(OBJ)
#	$(FC) $(FCFLAGS) -o $(EXE) $(OBJ)
#	@mv $(EXE) $(EXEDIR)
#	@ln -s $(EXEDIR)$(EXE) .

clean:
	@rm *.o

mod_smvgear_comode.o : mod_smvgear_comode.f
mod_smvgear_core.o : mod_smvgear_core.f mod_smvgear_comode.o
mod_smvgear_interface.o : mod_smvgear_interface.f mod_smvgear_comode.o mod_smvgear_core.o	
