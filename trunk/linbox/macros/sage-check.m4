# Copyright (c) the LinBox group
# This file is part of LinBox
# see COPYING for licence


AC_DEFUN([LB_CHECK_SAGE],
[
AC_MSG_CHECKING([whether to compile the SAGE interface])

AC_ARG_ENABLE(sage,
[AC_HELP_STRING([--enable-sage], [Enable the compilation of the SAGE interface])])
AM_CONDITIONAL(LINBOX_HAVE_SAGE, test "x$enable_sage" = "xyes")
AS_IF([test "x$enable_sage" = "xyes"],
[
AC_MSG_RESULT(yes)

if test "x$HAVE_NTL" = "xyes" ; then
	dnl  AC_CHECK_TOOL(OBJDUMP, objdump, false)
	AC_MSG_CHECKING([whether NTL was built with -fPIC])
	res=yes;
	if test -f "$NTL_HOME/lib/libntl.a" ; then
		$OBJDUMP --reloc $NTL_HOME/lib/libntl.a | $EGREP '(GOT|PLT|JU?MP_SLOT)' >/dev/null || res=no
	else if  test -f "$NTL_HOME/lib/libntl.so" ; then
			$OBJDUMP -R $NTL_HOME/lib/libntl.so | $EGREP '(GOT|PLT|JU?MP_SLOT)' >/dev/null || res=no
		else
			AC_MSG_RESULT(no, libntl not found !)
		fi
	fi

	if test "x$res" = "xno" ; then
		AC_MSG_RESULT(no)
		echo
		echo "You must have NTL compiled with -fPIC for Sage interface  "
		exit 1
	else
		AC_MSG_RESULT(yes)
	fi
fi

],[
AC_MSG_RESULT(no)
])
])
