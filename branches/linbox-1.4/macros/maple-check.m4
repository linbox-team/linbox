dnl Check for Maple, for maple examples
dnl Rich Seagraves, 2002-10-31
dnl Modified by PAscal Giorgi, 2003-12-03
dnl a gutted and modified version of ntl-check.m4 by Bradford Hovinen, 2001-06-13
dnl which was in turn inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
dnl Stolen from Chris Lahey
dnl stolen from Manish Singh again
dnl stolen back from Frank Belew
dnl stolen from Manish Singh
dnl Shamelessly stolen from Owen Taylor


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




#Enables the compilation of the drivers
AC_DEFUN([LB_DRIVER],
[
AC_MSG_CHECKING([whether to compile the drivers])

AC_ARG_ENABLE(drivers, [AC_HELP_STRING([--enable-drivers], [Enable the compilation of the drivers])],
[
AC_MSG_RESULT(yes)
compile_drivers="yes"
],[
AC_MSG_RESULT(no)
compile_drivers="no"
])
AM_CONDITIONAL(LINBOX_COMPILE_DRIVERS, test "x$compile_drivers" = "xyes")
])

dnl LB_CHECK_MAPLE()
AC_DEFUN([LB_CHECK_MAPLE],
[

AC_ARG_WITH(maple,
[AC_HELP_STRING([--with-maple=<path>|yes|no], [Use Maple library. If argument is no, you do not
  			    have the library installed on your machine (set as
			    default). If argument is yes or <empty> that means
			    the library is well installed and so reachable.
			    Otherwise you give the <path> to the directory which
			    contains the Software.
	    ])],
            [if test "$withval" = yes ; then
		MAPLE_HOME_PATH="${DEFAULT_CHECKING_PATH} unknown"
	      elif test "$withval" != no ; then
		MAPLE_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH} unknown"
	     fi	],
	    [])

AC_ARG_ENABLE(shared, dnl this name is a bad idea
[AC_HELP_STRING([--enable-shared], [Check for shared compilation (needed by --with-maple)])],
	      [have_shared="$enableval"],
              [have_shared="no"])

if test -n "$MAPLE_HOME_PATH" ; then
min_maple_version=ifelse([$1], ,9.0,$1)
AC_MSG_CHECKING(for MAPLE >= $min_maple_version)
fi

eval LIB_DIR=${libdir}
if test ${LIB_DIR} = "NONE/lib" ; then
	eval LIB_DIR="${prefix}/lib"
fi


for MAPLE_HOME in ${MAPLE_HOME_PATH}
do
		if test "x$MAPLE_HOME" != "xunknown"; then
			if test -r "${MAPLE_HOME}/bin/maple.system.type" && test -r "${MAPLE_HOME}/extern/include/maplec.h" ; then
				MAPLE_BIN=${MAPLE_HOME}/`${MAPLE_HOME}/bin/maple.system.type`
			else
				MAPLE_BIN=""
			fi
		else
			if test -r "/usr/local/bin/xmaple" && test -r "/usr/local/bin/maple.system.type"; then
					MAPLE_HOME=`sed -ne "s/MAPLE='\(.*\)'/\\1/p" /usr/local/bin/xmaple`
					MAPLE_BIN="${MAPLE_HOME}/"`${MAPLE_HOME}/bin/maple.system.type`
			elif test -r "/usr/bin/xmaple" ; then
					MAPLE_HOME=`sed -ne "s/MAPLE='\(.*\)'/\\1/p" /usr/bin/xmaple`
					MAPLE_BIN="${MAPLE_HOME}/"`${MAPLE_HOME}/bin/maple.system.type`
			else
					MAPLE_BIN=""
			fi
		fi

		if test -z "${MAPLE_BIN}" ; then
			maple_found="no"
		else
			maple_found="yes"
			if test $have_shared = "yes"; then
				${MAPLE_HOME}/bin/maple macros/maple-check-version.mpl > /dev/null
				MAPLE_VERSION=`cat maple_version.txt`

				if test ${MAPLE_VERSION} -lt 9; then
					AC_MSG_RESULT(problem)
					echo " your version of Maple is too old, at least version 9 is recquired. Disabling."
					break
				else
					AC_MSG_RESULT(found)
					LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${MAPLE_BIN}:${LIB_DIR}"
					LD_RUN_PATH="${LD_RUN_PATH}:${MAPLE_BIN}:${LIB_DIR}"
					export LD_LIBRARY_PATH
					export LD_RUN_PATH
					MAPLE_LIBS="-L${MAPLE_BIN} -lmaplec -lstdc++"
					MAPLE_CFLAGS="-I${MAPLE_HOME}/extern/include"
       					AC_SUBST(MAPLE_LIBS)
					AC_SUBST(MAPLE_CFLAGS)
					AC_SUBST(MAPLE_HOME)
					AC_SUBST(MAPLE_VERSION)
					if test ${MAPLE_VERSION} -ge 10 	; then
						AC_DEFINE_UNQUOTED(MAPLE_GMP_ACCESS,, [define is the version of Maple have access function to gmp data])
					fi
					AC_DEFINE(HAVE_MAPLE,1,[Define if MAPLE is installed])
					HAVE_MAPLE=yes
					break
				fi
			else
				AC_MSG_RESULT(problem)
				echo " you need to give option --enable-shared to allow Maple interfacing. Disabling."
				break
			fi
		fi
done

if test "x$maple_found" = "xno" ; then
	AC_MSG_RESULT(not found)
	AC_DEFINE(HAVE_MAPLE,0)
fi

AM_CONDITIONAL(LINBOX_HAVE_MAPLE, test "x$HAVE_MAPLE" = "xyes")
AM_CONDITIONAL(LINBOX_COMPILE_DRIVERS, test "x$compile_drivers" = "xyes" -o "x$HAVE_MAPLE" = "xyes" )

])

