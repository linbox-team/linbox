# Check for ALP
# Pascal Giorgi, 2001-10-10
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_ALP ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for the ALP library and define ALP_CFLAGS and ALP_LIBS

AC_DEFUN([LB_CHECK_ALP],
[

AC_ARG_WITH(alp-prefix,[  --with-alp-prefix=PFX      Prefix where ALP is installed (optional)],
[alp_prefix="$withval"],[alp_prefix=""])

min_alp_version=ifelse([$1], ,1.0 ,$1)
AC_MSG_CHECKING(for ALP >= $min_alp_version)



if test x$alp_prefix != x; then
	LD_LIBRARY_PATH=$alp_prefix/lib:$LD_LIBRARY_PATH
	export LD_LIBRARY_PATH
	CPLUS_INCLUDE_PATH=$alp_prefix/include:$CPLUS_INCLUDE_PATH
	export CPLUS_INCLUDE_PATH

	AC_SUBST(alp_prefix)
	ALP_CFLAGS="-I$(alp_prefix)/include"
	ALP_LIBS="-L$(alp_prefix)/lib -lalp"
	AC_SUBST(ALP_CFLAGS)
	AC_SUBST(ALP_LIBS)
	AC_DEFINE(HAVE_ALP)
	AC_MSG_RESULT(found)
fi

echo ALP IS AT ............. $alp_prefix
dnl Check for existence

AC_CHECK_LIB(alp, GetTime,
		  [echo alp check succeeded],[echo alp check for GetTime failed])
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
