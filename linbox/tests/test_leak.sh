#!/bin/bash - 

# Copyright (c) 2010 the LinBox group
# written by Brice Boyer <brice.boyer@imag.fr>
# This file is part of LinBox
# ========LICENCE========
# This file is part of the library FFLAS-FFPACK.
#
# FFLAS-FFPACK is free software: you can redistribute it and/or modify
# it under the terms of the  GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========
#/




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


