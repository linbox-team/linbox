include ../make.gen

#######################################################################
#
#	Gen_SuperLU package (version 1.0)
#	by Ahmet Duran and David Saunders
#	Computer & Information Science Department
#	University of Delaware
#	August 2002
#
#  To use the package:
#  1. Edit make.gen
#	to update LINBOXLIB and INCLUDESLINBOX definitions
#  2. Edit genlinsol.C to choose the Field such as
#	LinBox::UnparametricField<double>
#	LinBox::UnparametricField<float>
#	LinBox::Modular<LinBox::uint32>
#	LinBox::GivaroZpz<Std32>
#  3. Type
#	make
#  4. To run type, for example,
#	genlinsol < min5x5G
#  To remove the object files after the executable files have been
#  created, enter
#       make clean
#
#
#  We used SuperLU (Version 2.0) 
#  Thanks to
#  Univ. of California Berkeley, Xerox Palo Alto Research Center,
#  and Lawrence Berkeley National Lab.
#######################################################################

HEADER   = ../src

#SPIENV  = sp_ienv.c
GLINSOL = genlinsol.C
GRANK   = genrank.C 
GRANKN   = genrank_n.C
GRANKM   = genrank_m.C
GRANKG = genrank_g.C
GRANKG1 = genrank_g1.C
GRANKS   = genrank_simple.C
GENSUPERLU = gensuperlu.C

all:    field

field: genlinsol gensuperlu genrank genrank_n genrank_m genrank_simple \
	genrank_g genrank_g1

#sp_ienv.o: $(SPIENV) 
#	$(LOADER) $(LOADOPTS) $(SPIENV) \
#	-lm -c

genlinsol: $(GLINSOL) 
	$(LOADER) $(LOADOPTS) $(GLINSOL) \
	$(INCLUDESLINBOX) $(LINBOXLIB) \
	-lm -o genlinsol

gensuperlu: $(GENSUPERLU) 
	$(LOADER) $(LOADOPTS) $(GENSUPERLU) \
	$(INCLUDESLINBOX) $(LINBOXLIB) \
	-lm -o gensuperlu

genrank: $(GRANK) 
	$(LOADER) $(LOADOPTS) $(GRANK) \
	$(INCLUDESLINBOX) $(LINBOXLIB) \
	-lm -o genrank

genrank_n: $(GRANKN) 
	$(LOADER) $(LOADOPTS) $(GRANKN) \
	$(INCLUDESLINBOX) $(LINBOXLIB) \
	-lm -o genrank_n

genrank_m: $(GRANKM) 
	$(LOADER) $(LOADOPTS) $(GRANKM) \
	$(INCLUDESLINBOX) $(LINBOXLIB) \
	-lm -o genrank_m

genrank_g: $(GRANKG)
	$(LOADER) $(LOADOPTS) $(GRANKG) \
	$(INCLUDESLINBOX) $(LINBOXLIB) \
	-lm -o genrank_g

genrank_g1: $(GRANKG1)
	$(LOADER) $(LOADOPTS) $(GRANKG1) \
	$(INCLUDESLINBOX) $(LINBOXLIB) \
	-lm -o genrank_g1

genrank_simple: $(GRANKS) 
	$(LOADER) $(LOADOPTS) $(GRANKS) \
	$(INCLUDESLINBOX) $(LINBOXLIB) \
	-lm -o genrank_simple

F77EXM	= f77_main.o hbcode1.o c_bridge_dgssv.o

f77exm: $(F77EXM) ../$(SUPERLULIB)
	$(FORTRAN) $(F77EXM) ../$(SUPERLULIB) $(BLASLIB) -o $@

c_bridge_dgssv.o: c_bridge_dgssv.c
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

.c.o:
	$(CC) $(CFLAGS) -I$(HEADER) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

clean:	
	rm -f *.o *linsol *~ f77exm *superlu *rank *rank_m 
