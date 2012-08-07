dnl Check for OpenCL
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

dnl LB_CHECK_OCL([MINIMUM-VERSION [, ACTION_IF_FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for an OpenCL

AC_DEFUN([LB_CHECK_OCL],[
	AC_ARG_WITH(ocl,
	[AC_HELP_STRING([--with-ocl=<path>|yes],
	[Use OpenCL. This library is optional for LinBox compilation.
	 If argument is yes or <empty> that means the library is reachable
	 with the standard search path "/usr" or "/usr/local" (set as default).
	 Otherwise you give the <path> to the directory which contains the library.])],

	[if test "$withval" = yes; then
		OCL_HOME_PATH="/opt/AMDAPP ${DEFAULT_CHECKING_PATH}"
	elif test "$withval" != no; then
		OCL_HOME_PATH="$withval /opt/AMDAPP ${DEFAULT_CHECKING_PATH}"
	fi],

	[OCL_HOME_PATH="/opt/AMDAPP $DEFAULT_CHECKING_PATH}"])

	BACKUP_CXXFLAGS=${CXX_FLAGS}
	BACKUP_LIBS=${LIBS}

	AC_MSG_CHECKING(for OpenCL >= 1.0)

for OCL_HOME in ${OCL_HOME_PATH}
do
	if test -r "$OCL_HOME/include/CL/cl.h"; then
		if test "x$OCL_HOME" != "x/usr" -a "x$OCL_HOME" != "x/usr/local"; then
			OCL_CFLAGS="-I${OCL_HOME}/include"
			if test "x$OCL_HOME" != "x/opt/AMDAPP"; then
				OCL_LIBS="-L${OCL_HOME}/lib64 -L${OCL_HOME}/lib -lOpenCL -lpthread"
			else
				OCL_LIBS="-L${OCL_HOME}/lib/x86_64 -L${OCL_HOME}/lib/x86 -lOpenCL -lpthread"
			fi
		else
			OCL_CFLAGS=
			OCL_LIBS="-lOpenCL -lpthread"
		fi

		CXXFLAGS="${CXXFLAGS} ${OCL_CFLAGS}"
		LIBS="${LIBS} ${OCL_LIBS}"

		ocl_found="yes"
		break

	else
		ocl_found="no"
	fi
done

if test "x$ocl_found" = "xyes"; then
	AC_SUBST(OCL_CFLAGS)
	AC_SUBST(OCL_LIBS)
	AC_DEFINE(HAVE_OCL,1,[Define if OpenCL is installed])
	HAVE_OCL=yes

	AC_MSG_RESULT(found)
else
	AC_MSG_RESULT(not found)
fi

AM_CONDITIONAL(LINBOX_HAVE_OCL, test "x$HAVE_OCL" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
])
