# Copyright (c) the LinBox group
# This file is part of LinBox
# see COPYING for licence


AC_DEFUN([LB_CHECK_SAGE],
[
AC_MSG_CHECKING([whether to compile the sage interface])

AC_ARG_ENABLE(sage,
[AC_HELP_STRING([--enable-sage], [Enable the compilation of the sage interface])],
[
AC_MSG_RESULT(yes)
sage_interface="yes"
],[
AC_MSG_RESULT(no)
sage_interface="no"
])
AM_CONDITIONAL(LINBOX_HAVE_SAGE, test "x$sage_interface" = "xyes")
])
