include Makefile_Header
OUTPATH :=
EXE := smvgear.exe
EXEDIR := 
	
SRC := $(wildcard *.f)
OBJ := $(SRC:.f=.o)

all:$(OBJ)
	$(FC) $(FCFLAGS) -o $(EXE) $(OBJ)
#	@mv $(EXE) $(EXEDIR)
#	@ln -s $(EXEDIR)$(EXE) .

clean:
	@rm *.o

	
backsub.o   : backsub.f mod_comode.o
decomp.o    : decomp.f  mod_comode.o
jsparse.o   : jsparse.f mod_comode.o
ksparse.o   : ksparse.f mod_comode.o
mod_comode.o: mod_comode.f
pderiv.o    : pderiv.f  mod_comode.o
smvgear.o   : smvgear.f mod_comode.o
smvgear_init.o: smvgear_init.f mod_comode.o
smvgear_setchem.o : smvgear_setchem.f mod_comode.o
smvgear_setpara.o : smvgear_setpara.f mod_comode.o
smvgear_solve.o: smvgear_solve.f mod_comode.o
subfun.o    : subfun.f mod_comode.o
update.o    : update.f mod_comode.o
getdata.o : getdata.f
simpleintegrator.o : simpleintegrator.f
