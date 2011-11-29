AC_DEFUN([LB_CHECK_SAGE],
[
AC_MSG_CHECKING([whether to compile the sage interface])

AC_ARG_ENABLE(sage, 
[  --enable-sage Enable the compilation of the sage interface])

AM_CONDITIONAL(LINBOX_HAVE_SAGE, test "x$enable_sage" = "xyes")
AS_IF([test "x$enable_sage" = "xyes"],
[AC_MSG_RESULT(yes)],
[AC_MSG_RESULT(no)])
])

