#!/bin/bash - 

# Copyright (c) 2010 the LinBox group
# written by Brice Boyer <brice.boyer@imag.fr>
# This file is part of LinBox
# see COPYING for licence

# fast check (pwd = linbox/tests ; make check )
# for i in `find ./ -executable -type f | grep -v benchmark` ; do  valgrind $i ; done 2>&1 | less

# doc :
# recompiles a test-* and valgrinds it.

if [ $# -ne 1 ]
then
	echo "Usage: `basename $0` prog"
	exit -1
fi

if [ ! -f "${1}.C" ] 
then
	echo "not a test"
	exit -2
fi

if [ -f "${1}.o" ] 
then
	/bin/rm $1.o
else
	echo "pas de $1.o"
fi	

make $1 CXXFLAGS+="-O0 -g -DDEBUG" 

valgrind --leak-check=full --show-reachable=yes ./$1

exit 0

set -o nounset                              # Treat unset variables as an error


