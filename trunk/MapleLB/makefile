# LinBox/Maple interface makefile.  There are a few steps a user
# must take to install this library

# 1.  Please change ... to the directory in which you plan to install
# the interface.  Please be explicit, because this line will appear
# as a string literal inside the declaration of the Linbox module, and if it
# isn't written out, Maple will not be able to find the library and link to it
INSTALLDIR =... 

# 2.  As this is a package for use with Maple, the system will need to know
# where Maple is located on your system.  Please replace ... with the
# path of your maple installation
MAPLEPATH =... 

# 3.  Waterloo Maple has seen it wise to create a different directory for
# binaries in each platform.  Please replace bin.??? with the entry for
# your platform according to the following table:
# Solaris     - bin.SUN_SPARC_SOALRIS
# HP-UX       - bin.HP_RISC_UNIX
# IRIX        - bin.SGI_MIPS_UNIX
# AIX         - bin.IBM_RISC_UNIX
# OSF1/True64 - bin.DEC_ALPHA_UNIX
# Linux	      - bin.IBM_INTEL_LINUX
MAPLELIB =$(MAPLEPATH)/bin.??? 

# 4.  The interface depends on two other packages - LinBox and GMP.  Please
# specify the location of GMP's include directory.  Also, please specify
# the location of the LinBox distribution's main directory, ie the
# directory containing "linbox/blackbox" and "linbox/field" etc etc
GMPINCLUDE = ... 
LBPATH = ... 

# 5.  The compiler also needs to know the location of the libraries for the
# LinBox and GMP packages.  Please edit the ... in the following two lines
# to the appropriate entry
LBLIB = .../lib 
GMPLIB = .../lib 

# 6.  The compiler.  Please edit the following line to specify which C++
# compiler you would like to use.  I would highly recommend using
# g++ 3.0 or later (as my previous attempts to build the interface library
# with g++ 2.9 consistantly crashed Maple :-)
CC = g++

# 7.  Solaris users only.  If you are using a solaris box, you
#  will also need to specify the directory of your system's libstdc++ file.
#  If anybody knows a better way to do this, please email me, I beg you!
STDCPPLIB = .../lib/ 

## you shouldn't have to change anything below this point

CCFLAGS = -O2 -G

INCLUDES = -I. -I$(LBPATH) -I$(GMPINCLUDE) -I$(MAPLEPATH)/extern/include/

# Test for Solaris!  EW!
ifeq ($(MAPLELIB),$(MAPLEPATH)/bin.SUN_SPARC_SOLARIS)
	LIBS = -L$(LBLIB) -R$(LBLIB) -L$(GMPLIB) -R$(GMPLIB) -L$(MAPLELIB) \
               -R$(MAPLELIB) -R$(STDCPPLIB) -lmaplec -llinbox -lgmp 
else
	LIBS = -L$(LBLIB) -L$(GMPLIB) -L$(MAPLELIB) -lmaplec -llinbox -lgmp
endif

all: MapleLB.cc MapleBB.h 
	$(CC) $(CCFLAGS) MapleLB.cc -o libMapleLB.so -DDISABLE_COMMENTATOR $(INCLUDES) $(LIBS)

install: all linbox.mpl.head linbox.mpl.tail
	if test $(INSTALLDIR) != pwd ; then \
		mv libMapleLB.so $(INSTALLDIR)/libMapleLB.so; \
	fi; \
	cat linbox.mpl.head > $(INSTALLDIR)/linbox.mpl; \
	echo $(INSTALLDIR)"/libMapleLB.so\";\n" >> $(INSTALLDIR)/linbox.mpl; \
	cat linbox.mpl.tail >>$(INSTALLDIR)/linbox.mpl;

dist: clean 
	/usr/local/gnu/bin/tar czf ../MapleLB.tgz .

clean :
	rm -f libMapleLB.so $(INSTALLDIR)/libMapleLB.so linbox.mpl $(INSTALLDIR)/linbox.mpl *.tar *.gz *.tgz
