# Makefile for JX
#------------------------------------------------------------------------------

OBJ   = cmn_fjx70_mod.o fjx70_sub_mod.o fjx70_init_mod.o fjx70.o

FC    = ifort
FLAGS = -O2 -ip -fpp2 -W0 -assume byterecl

%.o: %.f
	$(FC) $(FLAGS) -c $<

jx    : $(OBJ)
	$(FC) $(FLAGS) -o JX70 $(OBJ)

#------------------------ Cleaning -----------------------------------

clean :
	/bin/rm fort.* *.o *.s *.mod core




PC compiler instructions (digital fortran)

> df -c  fjx_cmn_mod.f90

> df -c  fjx_sub_mod.f90

> df -c  fjx_init_mod.f90

              and if no other *.obj in directory:

> df     fjx_main.f90  *.obj
