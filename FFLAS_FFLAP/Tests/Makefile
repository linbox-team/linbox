#-------------------------------------------------------------------------
# To be completed with your own installation paths
BLASROOT          =  /usr/local/algebra/ATLAS
BLASARCH          =  SunOS_Stimpy
GIVAROROOT 	  =  /usr/local/algebra/givaro-dev
#-------------------------------------------------------------------------

OPT_ATLAS = -L$(BLASROOT)/lib/${BLASARCH}  -lcblas  -latlas #-lpthread


include ${GIVAROROOT}/bin/givaro-makefile 
CXXFLAGS += -g
CXXFLAGS += -w
CXXFLAGS += -I.. -I${BLASROOT}/include -I${BLASROOT}/include/${BLASARCH} -I${GIVAROROOT}/include	 


#for stimpy
OPTFLAGS+=-O3 -msupersparc
#for Sun workstations ( aphrodite, algol )
#OPTFLAGS+=-mcpu=ultrasparc -mtune=ultrasparc -fomit-frame-pointer -O3
#for vega
#OPTFLAGS+=-march=pentium4
OPTFLAGS+=-funroll-all-loops -frerun-loop-opt 
OPTFLAGS+=-fexpensive-optimizations -felide-constructors 
OPTFLAGS+=-fforce-addr 
OPTFLAGS+=-fforce-mem 
OPTFLAGS+=-fstrength-reduce
OPTFLAGS+=-fstrict-aliasing -ffast-math
OPT+=${OPTFLAGS}
CXX = /usr/local/gnu/gcc-3.2.1/bin/g++  ${OPT}

LDFLAGS += ${OPT_ATLAS} 
SRCDIR = ../FFLAP ../FFLAS ..

all: Test_CharPoly Test_LUdivine Test_Rank Test_fgemm dense_generator sparse_generator Legendre

Test_CharPoly: Test_CharPoly.C ../FFLAP.h ../FFLAS.h ../FFLAS/FFLAS_fgemm.inl ../FFLAP/FFLAP_LUdivine.inl ../FFLAP/FFLAP_flaswp.inl ../FFLAS/FFLAS_ftrsm.inl ../FFLAP/FFLAP_CharPoly.inl ../FFLAP/FFLAP_CharPoly_KGLU.inl ../FFLAP/FFLAP_MinPoly.inl ../FFLAP/FFLAP_MinPoly_construct.inl ../FFLAS/FFLAS_ftrsv.inl

Test_LUdivine: Test_LUdivine.C ../FFLAP.h ../FFLAS.h ../FFLAS/FFLAS_fgemm.inl ../FFLAP/FFLAP_LUdivine.inl ../FFLAP/FFLAP_flaswp.inl ../FFLAS/FFLAS_ftrsm.inl

Test_Rank: Test_Rank.C ../FFLAP.h ../FFLAS.h ../FFLAS/FFLAS_fgemm.inl ../FFLAP/FFLAP_LUdivine.inl ../FFLAP/FFLAP_flaswp.inl ../FFLAS/FFLAS_ftrsm.inl

Test_fgemm: ../FFLAS.h ../FFLAS/FFLAS_fgemm.inl Test_fgemm.C ../FFLAS/FFLAS_ftrsm.inl 

Legendre: Legendre.C ../FFLAP.h ../FFLAS.h ../FFLAS/FFLAS_fgemm.inl ../FFLAP/FFLAP_LUdivine.inl ../FFLAP/FFLAP_flaswp.inl ../FFLAS/FFLAS_ftrsm.inl

Legendre2: Legendre2.C ../FFLAP.h ../FFLAS.h ../FFLAS/FFLAS_fgemm.inl ../FFLAP/FFLAP_LUdivine.inl ../FFLAP/FFLAP_flaswp.inl ../FFLAS/FFLAS_ftrsm.inl

dense_generator:dense_generator.C

sparse_generator:sparse_generator.C

clean:
	rm -f *~ Test_CharPoly Test_LUdivine Test_Rank Test_fgemm dense_generator sparse_generator Legendre

