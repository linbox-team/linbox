
# To install, you have to specify C++ and C compilers,
# along with some compile-time options.
# This is done by setting the variables CXX, CC, CFLAGS, and CXXFLAGS
# below.  If you have a GNU compiler, the defaults should be OK.
# You can also change some settings in the file ../include/NTL/config.h 
# before starting the build---these settings affect the performance
# of NTL, not its correctness.


CXX = g++
# A C++ compiler, e.g., g++, CC, xlC

CC = gcc
# A C or C++ compiler, e. g., gcc, cc, xlc.
# There are a few components written in C which may be compiled
# under C++, but C compilers tend to generate slightly better code.
# The C++ and C compilers (if different) must generate
# compatible code.

CFLAGS =  -I../include -O2 -g -D_THREAD_SAFE -mv8
# Flags for the C compiler

CXXFLAGS = $(CFLAGS)
# Flags for the C++ compiler (usually the same)

# Flags for the linker
LDFLAGS = -lpthread

# Useful compile-time options (CFLAGS):
# -I../include -O2 -g -mv8 -qarch=ppc 

# Useful compile-time options (CXXFLAGS)
# -+

###############################################################

# -I../include
#   This is needed so that the compiler finds all the relevant include
#   files when building NTL.

 
# -O2 
#   This turns on the optimizer, which is essential for
#   reasonable performance.   It is not recommended to
#   use a higher level of optimization.
 
# -g (for debugging)
#   You need this if you want to run a debugger, like gdb.
#   Caution: with some compilers (but not GNU's), this will
#   turn off optimization (sometimes without warning).

# -mv8 (for non-ancient Sparc stations)
#   On SPARC-10s and later, make sure you include this flag.
#   This enables the use of fast integer multiplication.
#   On earlier SPARC stations, you should NOT use this flag,
#   but you should consider setting -DNTL_SINGLE_MUL (see below).

# -qarch=ppc (for xlc/C compilers on PowerPCs)
#   On PowerPCs, this is needed with the xlc/C compilers to get
#   access to the full instruction set.

# -+ (for xlc/C compilers)
#   Add this to the CXXFLAGS line when using the AIX compiler xlC.

#  That's it!  You can ignore everything else in this file!

##################################################################




# object files

O01=FFT.o FacVec.o GF2.o GF2E.o GF2EX.o GF2EXFactoring.o GF2X.o GF2X1.o
O02=$(O01) GF2XFactoring.o GF2XVec.o GetTime.o HNF.o IsFinite.o LLL.o LLL_FP.o
O03=$(O02) LLL_QP.o LLL_RR.o LLL_XD.o RR.o WordVector.o ZZ.o ZZVec.o
O04=$(O03) ZZX.o ZZX1.o ZZXCharPoly.o ZZXFactoring.o ZZ_p.o ZZ_pE.o ZZ_pEX.o
O05=$(O04) ZZ_pEXFactoring.o ZZ_pX.o ZZ_pX1.o ZZ_pXCharPoly.o ZZ_pXFactoring.o
O06=$(O05) fileio.o lip.o lzz_p.o lzz_pE.o lzz_pEX.o lzz_pEXFactoring.o
O07=$(O06) lzz_pX.o lzz_pX1.o lzz_pXCharPoly.o lzz_pXFactoring.o
O08=$(O07) mat_GF2.o mat_GF2E.o mat_RR.o mat_ZZ.o mat_ZZ_p.o
O09=$(O08) mat_ZZ_pE.o mat_lzz_p.o mat_lzz_pE.o mat_poly_ZZ.o
O10=$(O09) mat_poly_ZZ_p.o mat_poly_lzz_p.o pair_GF2EX_long.o
O11=$(O10) pair_GF2X_long.o pair_ZZX_long.o pair_ZZ_pEX_long.o
O12=$(O11) pair_ZZ_pX_long.o pair_lzz_pEX_long.o pair_lzz_pX_long.o
O13=$(O12) quad_float.o tools.o vec_GF2.o vec_GF2E.o vec_GF2XVec.o
O14=$(O13) vec_RR.o vec_ZZ.o vec_ZZVec.o vec_ZZ_p.o vec_ZZ_pE.o
O15=$(O14) vec_double.o vec_long.o vec_lzz_p.o vec_lzz_pE.o vec_quad_float.o
O16=$(O15) vec_vec_GF2.o vec_vec_GF2E.o vec_vec_RR.o vec_vec_ZZ.o
O17=$(O16) vec_vec_ZZ_p.o vec_vec_ZZ_pE.o vec_vec_long.o vec_vec_lzz_p.o
O18=$(O17) vec_vec_lzz_pE.o vec_xdouble.o xdouble.o
O19=$(O18) G_LLL_FP.o G_LLL_QP.o G_LLL_XD.o G_LLL_RR.o GenPrime.o

OBJ=$(O19)

# library source files


S01=FFT.c FacVec.c GF2.c GF2E.c GF2EX.c GF2EXFactoring.c GF2X.c GF2X1.c
S02=$(S01) GF2XFactoring.c GF2XVec.c HNF.c IsFinite.c LLL.c LLL_FP.c LLL_QP.c
S03=$(S02) LLL_RR.c LLL_XD.c RR.c WordVector.c ZZ.c ZZVec.c ZZX.c ZZX1.c
S04=$(S03) ZZXCharPoly.c ZZXFactoring.c ZZ_p.c ZZ_pE.c ZZ_pEX.c 
S05=$(S04) ZZ_pEXFactoring.c ZZ_pX.c ZZ_pX1.c ZZ_pXCharPoly.c
S06=$(S05) ZZ_pXFactoring.c fileio.c lip.c lzz_p.c lzz_pE.c lzz_pEX.c
S07=$(S06) lzz_pEXFactoring.c lzz_pX.c lzz_pX1.c
S08=$(S07) lzz_pXCharPoly.c lzz_pXFactoring.c mat_GF2.c mat_GF2E.c
S09=$(S08) mat_RR.c mat_ZZ.c mat_ZZ_p.c mat_ZZ_pE.c mat_lzz_p.c mat_lzz_pE.c
S10=$(S09) mat_poly_ZZ.c mat_poly_ZZ_p.c mat_poly_lzz_p.c pair_GF2EX_long.c
S11=$(S10) pair_GF2X_long.c pair_ZZX_long.c pair_ZZ_pEX_long.c
S12=$(S11) pair_ZZ_pX_long.c pair_lzz_pEX_long.c pair_lzz_pX_long.c
S13=$(S12) quad_float.c tools.c vec_GF2.c vec_GF2E.c vec_GF2XVec.c vec_RR.c
S14=$(S13) vec_ZZ.c vec_ZZVec.c vec_ZZ_p.c vec_ZZ_pE.c vec_double.c
S15=$(S14) vec_long.c vec_lzz_p.c vec_lzz_pE.c vec_quad_float.c
S16=$(S15) vec_vec_GF2.c vec_vec_GF2E.c vec_vec_RR.c vec_vec_ZZ.c
S17=$(S16) vec_vec_ZZ_p.c vec_vec_ZZ_pE.c vec_vec_long.c vec_vec_lzz_p.c
S18=$(S17) vec_vec_lzz_pE.c vec_xdouble.c xdouble.c
S19=$(S18) G_LLL_FP.c G_LLL_QP.c G_LLL_XD.c G_LLL_RR.c GenPrime.c

SRC = $(S19)


# library header files

IN01= FFT.h FacVec.h GF2.h GF2E.h GF2EX.h GF2EXFactoring.h GF2X.h
IN02=$(IN01) GF2XFactoring.h GF2XVec.h HNF.h IsFinite.h LLL.h PrimeSource.h
IN03=$(IN02) RR.h WordVector.h ZZ.h ZZVec.h ZZX.h ZZXFactoring.h ZZ_p.h
IN04=$(IN03) ZZ_pE.h ZZ_pEX.h ZZ_pEXFactoring.h ZZ_pX.h ZZ_pXFactoring.h
IN05=$(IN04) fileio.h lip.h lzz_p.h lzz_pE.h lzz_pEX.h lzz_pEXFactoring.h
IN06=$(IN05) lzz_pX.h lzz_pXFactoring.h mat_GF2.h mat_GF2E.h mat_RR.h
IN07=$(IN06) mat_ZZ.h mat_ZZ_p.h mat_ZZ_pE.h mat_lzz_p.h mat_lzz_pE.h
IN08=$(IN07) mat_poly_ZZ.h mat_poly_ZZ_p.h mat_poly_lzz_p.h matrix.h
IN09=$(IN08) pair.h vector.h pair_GF2EX_long.h pair_GF2X_long.h
IN10=$(IN09) pair_ZZX_long.h pair_ZZ_pEX_long.h pair_ZZ_pX_long.h
IN11=$(IN10) pair_lzz_pEX_long.h pair_lzz_pX_long.h quad_float.h
IN12=$(IN11) tools.h vec_GF2.h vec_GF2E.h vec_GF2XVec.h vec_RR.h
IN13=$(IN12) vec_ZZ.h vec_ZZVec.h vec_ZZ_p.h vec_ZZ_pE.h vec_double.h
IN14=$(IN13) vec_long.h vec_lzz_p.h vec_lzz_pE.h vec_quad_float.h
IN15=$(IN14) vec_vec_GF2.h vec_vec_GF2E.h vec_vec_RR.h vec_vec_ZZ.h
IN16=$(IN15) vec_vec_ZZ_p.h vec_vec_ZZ_pE.h vec_vec_long.h vec_vec_lzz_p.h
IN17=$(IN16) vec_vec_lzz_pE.h vec_xdouble.h xdouble.h config.h

INCL=$(IN17)



# test data

TD1=BerlekampTestIn BerlekampTestOut CanZassTestIn CanZassTestOut 
TD2=$(TD1) ZZXFacTestIn ZZXFacTestOut MoreFacTestIn LLLTestIn LLLTestOut RRTestIn RRTestOut 
TD3=$(TD2) MatrixTestIn MatrixTestOut CharPolyTestIn  
TD4=$(TD3) CharPolyTestOut QuadTestIn QuadTestOut

TD = $(TD4)

# test source files

TS1=Test.c BerlekampTest.c CanZassTest.c ZZXFacTest.c MoreFacTest.c LLLTest.c
TS2=$(TS1) subset.c MatrixTest.c CharPolyTest.c RRTest.c QuadTest.c
TS3=$(TS2) GF2EXTest.c BitMatTest.c ZZ_pEXTest.c lzz_pEXTest.c TestThreads.c

TS = $(TS3)

# scripts

SCRIPTS=MakeGetTime TestScript dosify unixify

# auxilliary source

MD=MakeDesc.c MakeDescAux.c newnames.c

GT=GetTime1.c GetTime2.c GetTime3.c GetTime4.c GetTime5.c TestGetTime.c


# documentation


D01=GF2.txt GF2E.txt GF2EX.txt GF2EXFactoring.txt GF2X.txt
D02=$(D01) GF2XFactoring.txt GF2XVec.txt HNF.txt LLL.txt RR.txt
D03=$(D02) ZZ.txt ZZVec.txt ZZX.txt ZZXFactoring.txt ZZ_p.txt ZZ_pE.txt
D04=$(D03) ZZ_pEX.txt ZZ_pEXFactoring.txt ZZ_pX.txt ZZ_pXFactoring.txt
D05=$(D04) conversions.txt flags.txt lzz_p.txt lzz_pE.txt lzz_pEX.txt
D06=$(D05) lzz_pEXFactoring.txt lzz_pX.txt lzz_pXFactoring.txt mat_GF2.txt
D07=$(D06) mat_GF2E.txt mat_RR.txt mat_ZZ.txt mat_ZZ_p.txt mat_ZZ_pE.txt
D08=$(D07) mat_lzz_p.txt mat_lzz_pE.txt mat_poly_ZZ.txt mat_poly_ZZ_p.txt
D09=$(D08) mat_poly_lzz_p.txt matrix.txt pair.txt vector.txt
D10=$(D09) quad_float.txt sedscript.txt tools.txt vec_GF2.txt
D11=$(D10) vec_GF2E.txt vec_RR.txt vec_ZZ.txt vec_ZZ_p.txt vec_ZZ_pE.txt
D12=$(D11) vec_lzz_p.txt vec_lzz_pE.txt xdouble.txt names.txt
D13=$(D12) tour-ack.html tour-intro.html tour-time.html tour-changes.html
D14=$(D13) tour-modules.html tour-unix.html tour-examples.html
D15=$(D14) tour-roadmap.html tour-win.html tour-impl.html tour-struct.html
D16=$(D15) tour.html

DOC = $(D16)



# test program executables

PROG1=Test BerlekampTest CanZassTest ZZXFacTest MoreFacTest LLLTest  BitMatTest
PROG2=$(PROG1) MatrixTest CharPolyTest RRTest QuadTest 
PROG3=$(PROG2) GF2EXTest subset ZZ_pEXTest lzz_pEXTest TestThreads

PROGS = $(PROG3)

# saved to TAR file

SAVE1=README  makefile $(SRC) $(INCL) $(DOC)
SAVE2=$(SAVE1) $(SCRIPTS) $(MD) $(GT) $(TS) $(TD) mach_desc.win
SAVE = $(SAVE2)

SFILES=makefile $(SRC) $(SCRIPTS) $(MD) $(GT) $(TS) $(TD) mach_desc.win

DIRNAME=ntl-3.7a
WINDIR=WinNTL-3_7a

install:
	$(CC) -c $(CFLAGS) MakeDescAux.c
	$(CC) -o MakeDesc $(CFLAGS) MakeDesc.c MakeDescAux.o
	./MakeDesc
	mv mach_desc.h ../include/NTL/mach_desc.h
	make ntl.a
	make Test
	./Test
	-rm Test
	sh TestScript


ntl.tar: 
	tar -cvf ntl.tar $(SAVE)

package:
	sh unixify "$(SFILES)" "$(INCL)" "$(DOC)" 
	-rm -r $(DIRNAME)
	mv unix $(DIRNAME)
	tar -cvf $(DIRNAME).tar $(DIRNAME)
	gzip $(DIRNAME).tar
	- rm -r $(DIRNAME)

winpack:
	sh dosify "$(SRC)" "$(INCL)" "$(DOC)" "$(TS)" "$(TD)"
	- rm -r $(WINDIR)
	mv dos $(WINDIR)
	find ./$(WINDIR) -print | zip -l $(WINDIR) -@
	rm -r $(WINDIR)

ntl.a:	$(OBJ) 
	ar ruv ntl.a $(OBJ) /usr/lib/libefence.a
	- ranlib ntl.a


clobber:	
	- rm ntl.a mach_desc.h ../include/NTL/mach_desc.h MakeDesc GetTime.c TestGetTime 
	- rm $(PROGS)
	- rm *.o

clean:
	- rm *.o 
	- rm MakeDesc TestGetTime

lip.o:	lip.c 
	$(CC) $(CFLAGS) -c lip.c

IsFinite.o:	IsFinite.c
	$(CC) $(CFLAGS) -c IsFinite.c


GetTime.c:
	sh MakeGetTime $(CC) $(CFLAGS)

GetTime.o: GetTime.c 
	$(CC) $(CFLAGS) -c GetTime.c


.c.o: 
	$(CXX) $(CXXFLAGS) -c $<

.c: 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< ntl.a  -lm

