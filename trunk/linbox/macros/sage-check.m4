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


AC_DEFUN([LB_CHECK_SAGE],
[
AC_MSG_CHECKING([whether to compile the SAGE interface])

AC_ARG_ENABLE(sage,
[AC_HELP_STRING([--enable-sage], [Enable the compilation of the SAGE interface])])
AM_CONDITIONAL(LINBOX_HAVE_SAGE, test "x$enable_sage" = "xyes")
AS_IF([test "x$enable_sage" = "xyes"],
[
AC_MSG_RESULT(yes)

dnl  if test "x$HAVE_NTL" = "xyes" ; then
	dnl  AC_CHECK_TOOL(OBJDUMP, objdump, false)
	dnl  AC_MSG_CHECKING([whether NTL was built with -fPIC])
	dnl  res=yes;
	dnl  if test -f "$NTL_HOME/lib/libntl.a" ; then
		dnl  $OBJDUMP --reloc $NTL_HOME/lib/libntl.a | $EGREP '(GOT|PLT|JU?MP_SLOT)' >/dev/null || res=no
	dnl  else if  test -f "$NTL_HOME/lib/libntl.so" ; then
			dnl  $OBJDUMP -R $NTL_HOME/lib/libntl.so | $EGREP '(GOT|PLT|JU?MP_SLOT)' >/dev/null || res=no
		dnl  else
			dnl  AC_MSG_RESULT(no, libntl not found !)
		dnl  fi
	dnl  fi

	dnl  if test "x$res" = "xno" ; then
		dnl  AC_MSG_RESULT(no)
		dnl  echo
		dnl  echo "You must have NTL compiled with -fPIC for Sage interface  "
	dnl	dnl  exit 1
		dnl  else
		dnl  AC_MSG_RESULT(yes)
	dnl  fi
dnl  fi

],[
AC_MSG_RESULT(no)
])
])
