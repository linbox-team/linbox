dnl Copyright (c) 2011 the LinBox group
dnl This file is part of LinBox
dnl see COPYING for licence


dnl aclocal-include.m4
dnl This macro adds the name macrodir to the set of directories
dnl that `aclocal' searches for macros.

dnl serial 1

dnl AM_ACLOCAL_INCLUDE(macrodir)
AC_DEFUN([AM_ACLOCAL_INCLUDE],
[
	AM_CONDITIONAL(INSIDE_GNOME_COMMON, test x = y)

	test -n "$ACLOCAL_FLAGS" && ACLOCAL="$ACLOCAL $ACLOCAL_FLAGS"

	for k in $1 ; do ACLOCAL="$ACLOCAL -I $k" ; done
])
