# Check for LIDIA
# Pascal Giorgi, 2001-12-10
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_LIDIA ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for the LIDIA library and define LIDIA_CFLAGS and LIDIA_LIBS

AC_DEFUN([LB_CHECK_LIDIA],
[

AC_ARG_WITH(lidia-prefix,[  --with-lidia-prefix=PFX      Prefix where LIDIA is installed (optional)],
[lidia_prefix="$withval"],[lidia_prefix=""])

min_lidia_version=ifelse([$1], ,2.1 ,$1)
AC_MSG_CHECKING(for LIDIA >= $min_lidia_version)



if test x$lidia_prefix != x; then
	LD_LIBRARY_PATH=$lidia_prefix/lib:$LD_LIBRARY_PATH
	export LD_LIBRARY_PATH
	CPLUS_INCLUDE_PATH=$lidia_prefix/include:$CPLUS_INCLUDE_PATH
	export CPLUS_INCLUDE_PATH

	AC_SUBST(lidia_prefix)
	LIDIA_CFLAGS="-I$(lidia_prefix)/include"
	LIDIA_LIBS="-L$(lidia_prefix)/lib -lLiDIA"
	AC_SUBST(LIDIA_CFLAGS)
	AC_SUBST(LIDIA_LIBS)
	AC_DEFINE(HAVE_LIDIA)
	AC_MSG_RESULT(found)
fi

echo LIDIA IS AT ............. $lidia_prefix
dnl Check for existence

AC_CHECK_LIB(lidia, GetTime,
		  [echo lidia check succeeded],[echo lidia check for GetTime failed])
#[
#dnl Check if the version is new enough
#dnl FIXME
#
#NTL_CFLAGS="-I$(ntl_prefix)/include"
#dnl NTL_LIBS="-L$(ntl_prefix)/src -lntl"
#NTL_LIBS="-L$(prefix)/lib -lntl"
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
