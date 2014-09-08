include ../Makefile_Header
	
SRC := $(wildcard *.f)
OBJ := $(SRC:.f=.o)

all:$(OBJ)

clean:
	@rm *.o

module_model_cheminfo.o : module_model_cheminfo.f module_model_parameter.o
module_model_parameter.o : module_model_parameter.f
