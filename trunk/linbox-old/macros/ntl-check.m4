# Check for NTL
# Bradford Hovinen, 2001-06-13
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_NTL ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for Victor Shoup's NTL (Number Theory Library) and define
dnl NTL_CFLAGS and NTL_LIBS

AC_DEFUN([LB_CHECK_NTL],
[

AC_ARG_WITH(ntl-prefix,[  --with-ntl-prefix=PFX      Prefix where NTL is installed (optional)],
[ntl_prefix="$withval"],[ntl_prefix=""])

min_ntl_version=ifelse([$1], ,4.0,$1)
AC_MSG_CHECKING(for NTL >= $min_ntl_version)

set echo
if test x$ntl_prefix != x; then
#	LD_LIBRARY_PATH=$ntl_prefix/src:$LD_LIBRARY_PATH
#	export LD_LIBRARY_PATH
# because NTL's library archive is not given name libntl.a, we try this:
	if test \( \! -d $prefix/lib \) ; then mkdir $prefix/lib ; fi
	ln -sn $ntl_prefix/src/ntl.a $prefix/lib/libntl.a
	LD_LIBRARY_PATH=$prefix/lib:$LD_LIBRARY_PATH
	export LD_LIBRARY_PATH
	CPLUS_INCLUDE_PATH=$ntl_prefix/include:$CPLUS_INCLUDE_PATH
	export CPLUS_INCLUDE_PATH
fi

dnl Check for existence

dnl FIXME -  Got tired of trying to make this check work, so skip it. -bds
#AC_CHECK_LIB(ntl, _ntl_GetTime,
#[
#dnl Check if the version is new enough
#dnl FIXME
#
#NTL_CFLAGS="-I$(ntl_prefix)/include"
#NTL_LIBS="-L$(ntl_prefix)/src -lntl"
#AC_SUBST(NTL_CFLAGS)
#AC_SUBST(NTL_LIBS)
#AC_DEFINE(HAVE_NTL)
#AC_MSG_RESULT(found)
#ifelse([$2], , :, [$2])
#],
#[
#AC_MSG_RESULT(not found)
#ifelse([$3], , :, [$3])
#])

])
