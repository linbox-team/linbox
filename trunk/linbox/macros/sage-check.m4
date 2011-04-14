# Copyright (c) the LinBox group
# This file is part of LinBox
# see COPYING for licence


AC_DEFUN([LB_CHECK_SAGE],
[
AC_MSG_CHECKING([whether to compile the sage interface])

AC_ARG_ENABLE(sage,
[  --enable-sage Enable the compilation of the sage interface],
[
AC_MSG_RESULT(yes)
sage_interface="yes"

if test "x$HAVE_NTL" = "xyes" ; then
dnl  AC_CHECK_TOOL(OBJDUMP, objdump, false)
AC_MSG_CHECKING([whether NTL was built with -fPIC])
res=yes;
$OBJDUMP --reloc $NTL_HOME/lib/libntl.a | $EGREP '(GOT|PLT|JU?MP_SLOT)' >/dev/null || res=no
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
sage_interface="no"
])
AM_CONDITIONAL(LINBOX_HAVE_SAGE, test "x$sage_interface" = "xyes")
])
