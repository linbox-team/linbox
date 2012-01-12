dnl Copyright (c) the LinBox group
dnl This file is part of LinBox
dnl see COPYING for licence


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
