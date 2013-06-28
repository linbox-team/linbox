dnl Check for GMP
dnl Copyright (c) the LinBox group
dnl This file is part of LinBox

 dnl ========LICENCE========
 dnl This file is part of the library LinBox.
 dnl
 dnl LinBox is free software: you can redistribute it and/or modify
 dnl it under the terms of the  GNU Lesser General Public
 dnl License as published by the Free Software Foundation; either
 dnl version 2.1 of the License, or (at your option) any later version.
 dnl
 dnl This library is distributed in the hope that it will be useful,
 dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
 dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 dnl Lesser General Public License for more details.
 dnl
 dnl You should have received a copy of the GNU Lesser General Public
 dnl License along with this library; if not, write to the Free Software
 dnl Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 dnl ========LICENCE========
 dnl
dnl Modified by Pascal Giorgi, 2003-12-03

dnl LB_CHECK_GMP ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for the GNU Multiprecision library and define GMP_CFLAGS and GMP_LIBS

AC_DEFUN([LB_CHECK_GMP], [

		AC_ARG_WITH(gmp,
			[AC_HELP_STRING([--with-gmp=<path>|yes|no],[
				Use GMP library.
				If argument is no, you do not have the library installed on your machine.
				If argument is yes or <empty> that means the library is reachable with the standard
				search path "/usr" or "/usr/local"  (set as default).
				Otherwise you give the <path> to the directory which contain the library.
				])],
		[if test "$withval" = yes ; then
			GMP_HOME_PATH="/usr"
	         elif test "$withval" != no ; then
			GMP_HOME_PATH="$withval"
	        fi],
			[GMP_HOME_PATH=" /usr"])

		min_gmp_version=ifelse([$1], ,4.0.0,$1)

		dnl Check for existence

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for GMP >= $min_gmp_version with cxx enabled)

GMP_PATH=
for GMP_HOME in ${GMP_HOME_PATH}
do
		if test "x$GMP_HOME" != "x/usr" -a "x$GMP_HOME" != "x/usr/local"; then
		if test -r "$GMP_HOME/include/gmp.h" ; then
			GMP_CFLAGS="-I${GMP_HOME}/include"
			GMP_PATH="-L${GMP_HOME}/lib"
			GMP_LIBS="-L${GMP_HOME}/lib -lgmp"
		elif test -r "$GMP_HOME/gmp.h" ; then
			GMP_CFLAGS="-I${GMP_HOME}"
			GMP_PATH="-L${GMP_HOME}"
			GMP_LIBS="-L${GMP_HOME} -lgmp"
		else
			echo "($GMP_HOME) seems an invalid GMP prefix"
			echo "Searching GMP in PATH"
			GMP_CFLAGS=""
			GMP_LIBS="-lgmp"
		fi
	else
		GMP_CFLAGS=""
		GMP_LIBS="-lgmp"
		fi

		CXXFLAGS="${CXXFLAGS} ${GMP_CFLAGS}"
		LIBS="${LIBS} ${GMP_LIBS}"

		AC_TRY_LINK(
		[
		#ifdef __PATHCC__
		#include <rw/_defs.h>
		#include <ansi/_cstddef.h>
		#endif
		#include "stddef.h"
		#include <gmp.h>
		],
		[mpz_t a; mpz_init (a);],
		[
		dnl  # See if we are running GMP 4.0 with --enable-cxx
        		AC_TRY_RUN(
			[#include <gmpxx.h>
			int main () {
			if (__GNU_MP_VERSION < 4) return -1;
			mpz_class a(2),b(3),c(5); if ( a+b == c ) return 0; else return -1; }
			],
 			[
				AC_MSG_RESULT(found)
			AC_DEFINE(HAVE_GMP,1,[Define if GMP is installed])

			dnl  GMP_VERSION="" dnl I could find it but why is it here ?
			GMP_LIBS="${GMP_PATH} -lgmpxx -lgmp"
			dnl  AC_SUBST(GMP_VERSION)
			AC_SUBST(GMP_LIBS)
				AC_SUBST(GMP_CFLAGS)
			break;
	   			],[
			AC_MSG_RESULT(no : gmp is too old or not compiled with cpp support)
			dnl  AC_DEFINE(GMP_NO_CXX,1,[Define if GMP has no <gmpxx.h>])
			dnl  GMP_VERSION="-DGMP_NO_CXX"
			dnl  AC_SUBST(GMP_VERSION)
			],[ dnl This should never happen
					AC_MSG_RESULT(no)
				])
			],[
				AC_MSG_RESULT(unknown)
				echo "WARNING: You appear to be cross compiling, so there is no way to determine"
				echo "whether your GMP version is new enough. I am assuming it is."
				AC_SUBST(GMP_CFLAGS)
				AC_SUBST(GMP_LIBS)
				AC_DEFINE(HAVE_GMP,1,[Define if GMP is installed])
			])
		unset GMP_CFLAGS
		unset GMP_LIBS
done

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
