# Copyright(c)'2011 by The Givaro group
# Written by BB <bboyer@imag.fr>
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

dnl enable basic debug mode.
AC_DEFUN([AC_DEBUG],
[AC_MSG_CHECKING([whether to enable debugging options in the library])
  AC_ARG_ENABLE(debug,
[AC_HELP_STRING([--enable-debug], [enable debugging options in library])],
      USE_DEBUG=$enableval,
      USE_DEBUG=no)
  AC_MSG_RESULT([$USE_DEBUG])
  AM_CONDITIONAL(DEBUG, [test $USE_DEBUG = yes])
  DBG=$USE_DEBUG
  AC_SUBST(DBG)dnl
]
)

AC_DEFUN([AC_PROFILE],
[AC_MSG_CHECKING([whether to enable profiling everything in the library])
  AC_ARG_ENABLE(profile,
[AC_HELP_STRING([--enable-profile], [enable profiling options in library])],
      USE_PROFILE=$enableval,
      USE_PROFILE=no)
  AC_MSG_RESULT([$USE_PROFILE])
  AM_CONDITIONAL(PROFILE, [test $USE_PROFILE = yes])
  PROF=$USE_PROFILE
  AC_SUBST(PROF)dnl
]
)

dnl Enable warnings from compiler.
AC_DEFUN([AC_WARNINGS],
[AC_MSG_CHECKING([whether to enable warnings when compiling the library])
  AC_ARG_ENABLE(warnings,
[AC_HELP_STRING([--enable-warnings],  [enable warings when compiling the library])],
      USE_WARNINGS=$enableval,
      USE_WARNINGS=no)
  AC_MSG_RESULT([$USE_WARNINGS])
  AM_CONDITIONAL(WARNINGS, [test $USE_WARNINGS = yes])
  WARN=$USE_WARNINGS
  AC_SUBST(WARN)dnl
]
)

AC_DEFUN([AC_COMPILER_NAME],
[
AC_MSG_CHECKING(for family name of compiler)
AC_TRY_RUN(dnl ICC ?
[   #ifdef __INTEL_COMPILER
   int main() { return 0 ; }
   #else
   pas intel
   #endif],dnl
[dnl
   AC_MSG_RESULT(icc)
   CCNAM=icc
   AC_SUBST(CCNAM)
],dnl GCC ?
[dnl
   AC_TRY_RUN(dnl GCC ?
[#ifdef __GNUC__
   int main() { return 0 ; }
   #else
   pas gcc non plus.
   #endif],[
   AC_MSG_RESULT(gcc)
   CCNAM=gcc
   AC_SUBST(CCNAM)
   ],[
   AC_MSG_RESULT(unknown)
   CCNAM=unknown
   AC_SUBST(CCNAM)
   ],[
   AC_MSG_RESULT(unknown)
   CCNAM=unknown
   AC_SUBST(CCNAM)
   ])
],dnl GCC !
[
   AC_MSG_RESULT(unknown)
   CCNAM=unknown
   AC_SUBST(CCNAM)
])
])

