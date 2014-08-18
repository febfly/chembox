ROOTDIR = ..
LIB     = $(ROOTDIR)/lib
MOD     = $(ROOTDIR)/mod
RUN     = ../../run
G2R     = chemrxn_g2r.exe
G2K     = chemrxn_g2k.exe
RCA     = chemrxn_rca.exe
FPY     = chemrxn_fpy.exe
COT     = chemrxn_cot.exe

include $(ROOTDIR)/Makefile_header

SRC = module_photolysis.f module_spec.f module_reaction.f module_geos.f module_ream.f module_kpp.f module_analyze.f module_forpython.f
SE1 = dr_geos2ream.f
SE2 = dr_geos2kpp.f
SE3 = dr_ream_cascade.f
SE4 = dr_forpython.f
SE5 = dr_chemout.f

OBJ = $(SRC:.f=.o)
OE1 = $(SE1:.f=.o)
OE2 = $(SE2:.f=.o)
OE3 = $(SE3:.f=.o)
OE4 = $(SE4:.f=.o)
OE5 = $(SE5:.f=.o)


all : $(OBJ) $(OE1) $(OE2) $(OE3) $(OE4) $(OE5)
	$(AR) crs libchemrxn.a $(OBJ)
	mv libchemrxn.a $(LIB)
	$(COMPILER) $(LDFLAGS) $(OBJ) $(OE1) $(LINK) -o $(G2R)
	$(COMPILER) $(LDFLAGS) $(OBJ) $(OE2) $(LINK) -o $(G2K)
	$(COMPILER) $(LDFLAGS) $(OBJ) $(OE3) $(LINK) -o $(RCA)
	$(COMPILER) $(LDFLAGS) $(OBJ) $(OE4) $(LINK) -o $(FPY)
	$(COMPILER) $(LDFLAGS) $(OBJ) $(OE5) $(LINK) -o $(COT)
	cp -rf $(G2R) $(RUN)
	cp -rf $(G2K) $(RUN)
	cp -rf $(RCA) $(RUN)
	cp -rf $(FPY) $(RUN)
#	mv *.mod $(MOD)
clean :
	rm -f *.o *.mod

