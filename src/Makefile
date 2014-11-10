include Makefile_Header
OUTPATH :=
	
SRC := $(wildcard *.f)
OBJ := $(SRC:.f=.o)

all:$(OBJ)
	@$(MAKE) -C model_util all
	@$(MAKE) -C model_common all
	@$(MAKE) -C chemmech_geoschem all 
	@$(MAKE) -C chemmech_ream all 
	@$(MAKE) -C chemmech_mcm all
	@$(MAKE) -C solver_smvgear all 
	@$(MAKE) -C photo_fastjx all 
	@$(MAKE) -C model_main all 
#	$(FC) $(FCFLAGS) -o $(EXE) $(OBJ)
#	@mv $(EXE) $(EXEDIR)
#	@ln -s $(EXEDIR)$(EXE) .

clean:
	@$(MAKE) -C model_util clean
	@$(MAKE) -C model_common clean
	@$(MAKE) -C chemmech_geoschem clean
	@$(MAKE) -C chemmech_ream clean
	@$(MAKE) -C solver_smvgear clean 
	@$(MAKE) -C photo_fastjx clean 
	@$(MAKE) -C model_main clean 
	
