CC = g++

LBINCLUDE = 
GMPINCLUDE = 
MAPLEINCLUDE =

LBLIB =
GMPLIB =
MAPLELIB = 



.cc .so:
	/usr/local/gnu/gcc3/bin/g++ -G $*.cc -o lib$*.so -DDISABLE_COMMENTATOR -I/m/ecua/seagrave/linbox/linbox-new -I/usr/local/algebra/gmp/include/ -I/usr/local/maple/extern/include/ -R/usr/local/algebra/gmp/lib/ -R/usr/local/maple/bin.SUN_SPARC_SOLARIS/ -lmaplec -L/m/ecua/seagrave/linbox/linbox-new/lib/ -R/m/ecua/seagrave/linbox/linbox-new/lib/ -llinbox -lgmp -R/usr/local/gnu/gcc3/lib/

#%:%.cc 
#	g++ $@.cc -Ilinbox-new -I/usr/local/algebra/gmp/include
#-Llinbox-new/lib -Rlinbox-new/lib -L/usr/local/algebra/gmp/lib
#-R/usr/local/algebra/gmp/lib -R/usr/local/gnu/gcc/lib

server:
	/usr/local/gnu/gcc3/bin/g++ server.cc -o LBserver -DDISABLE_COMMENTATOR -I/m/ecua/seagrave/linbox/linbox-new/ -I/usr/local/algebra/gmp/include/ -L/usr/local/algebra/gmp/lib -R/usr/local/algebra/gmp/lib -L/m/ecua/seagrave/linbox/linbox-new/lib/ -R/m/ecua/seagrave/linbox/linbox-new/lib/ -llinbox -lgmp

clean :
	rm *.oo
