# Check for expat
# Rich Seagraves
# stolen from Pascal Giorgi, 2001-12-10
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_EXPAT ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for expat Library and define EXPAT_CFLAGS and EXPAT_LIBS

AC_DEFUN([LB_CHECK_EXPAT],
[

AC_ARG_WITH(expat-prefix,[  --with-expat-prefix=PFX Prefix where expat is installed (optional)],
[expat_prefix="$withval"],[expat_prefix=""])

min_expat_version=ifelse([$1], ,1.95,$1)
AC_MSG_CHECKING(for expat >= $min_expat_version)

if test "x${expat_prefix}" = "x"; then
	AC_MSG_RESULT(not found)
	BACKUP_CXXFLAGS=${CXXFLAGS}
	BACKUP_LIBS=${LIBS}
	HAVE_EXPAT=no
else

dnl Check for existence

if test "x${expat_prefix}" != "x/usr" -a "x${expat_prefix}" != "x/usr/local"; then
	EXPAT_CFLAGS="-I${expat_prefix}/include"
	EXPAT_LIBS="-L${expat_prefix}/lib -lexpat"
else
	EXPAT_CFLAGS=
	EXPAT_LIBS=-lexpat
fi

# By default, these should be empty. We set them to include real data
# only if LIDIA is actually found.

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

CXXFLAGS="${CXXFLAGS} ${EXPAT_CFLAGS} ${GMP_CFLAGS}"
LIBS="${CXXFLAGS} ${EXPAT_LIBS} ${GMP_LIBS}"

AC_TRY_RUN(
[#include <expat.h>
int main () { 
  if(XML_MAJOR_VERSION < 1 
    || (XML_MAJOR_VERSION == 1 && XML_MINOR_VERSION < 95)) 
 	return -1;
  else
	return 0;
}
],[
AC_MSG_RESULT(found)
AC_SUBST(EXPAT_CFLAGS)
AC_SUBST(EXPAT_LIBS)

AC_DEFINE(XMLENABLED,1,[Define if Expat is installed])

# expat was found, so make sure tests and headers get included.

HAVE_EXPAT=yes

],[
AC_MSG_RESULT(not found)
dnl echo "Sorry, Expat wasn't found or too old.  Disabling XML reading & writing."
HAVE_EXPAT=no

unset EXPAT_CFLAGS
unset EXPAT_LIBS

],[
AC_MSG_RESULT(unknown)
echo "WARNING: You appear to be cross compiling, so there is no way to determine"
echo "whether your expat version is right. I am assuming it isn't."

HAVE_EXPAT=no

unset EXPAT_CFLAGS
unset EXPAT_LIBS

])

fi
AM_CONDITIONAL(HAVE_EXPAT, test "x$HAVE_EXPAT" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}

])
