# Check for SACLIB
# Bradford Hovinen, 2001-06-13
# Modified by Pascal Giorgi, 2003-12-03 
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_SACLIB ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for the GNU Multiprecision library and define SACLIB_CFLAGS and SACLIB_LIBS

AC_DEFUN([LB_CHECK_SACLIB],
[

AC_ARG_WITH(saclib,
[  --with-saclib=<path>|yes|no Use Saclib library. If argument is no, you do 
                              not have the library installed on your machine 
			      (set as default). If argument is yes or <empty> 
			      that means the library is reachable with the 
			      standard search path (/usr or /usr/local).
			      Otherwise you give the <path> to the directory 
			      which contain the library. 
	     ],
	     [if test "$withval" = yes ; then
			SACLIB_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	      elif test "$withval" != no ; then
			SACLIB_HOME_PATH="$withval"
	     fi],
	     [])

min_saclib_version=ifelse([$1], ,2.0,$1)


dnl Check for existence

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

if test -n  "$SACLIB_HOME_PATH"; then
AC_MSG_CHECKING(for SACLIB >= $min_saclib_version)
fi

for SACLIB_HOME in ${SACLIB_HOME_PATH} 
 do
if test -r "$SACLIB_HOME/include/saclib.h"; then

	if test "x$SACLIB_HOME" != "x/usr" -a "x$SACLIB_HOME" != "x/usr/local"; then
		SACLIB_CFLAGS="-I${SACLIB_HOME}/include"
		SACLIB_LIBS="-L${SACLIB_HOME}/lib -lsaclib"
	else
		SACLIB_CFLAGS=
		SACLIB_LIBS="-lsaclib"		
	fi
	
	CXXFLAGS="${BACKUP_CXXFLAGS} ${SACLIB_CFLAGS} ${GMP_CFLAGS}" 	LIBS="${BACKUP_LIBS} ${SACLIB_LIBS} ${GMP_LIBS}"

	AC_TRY_LINK(
	[#include <saclib.h>],
	[BDigit a;],
	[
	AC_TRY_RUN(
	[#include <saclib.h>
	int main () {  if ( __GNU_MP_VERSION < 3) return -1; else return 0; }
	],[
	saclib_found="yes"		
	break
	],[
	saclib_problem="$problem $SACLIB_HOME"	
	unset SACLIB_CFLAGS
	unset SACLIB_LIBS

	
	],[
	saclib_found="yes"
	saclib_cross="yes"
	break
	])
	],
	[
	saclib_found="no"
	saclib_checked="$saclib_checked $SACLIB_HOME"	
	unset SACLIB_CFLAGS
	unset SACLIB_LIBS
	])
else
	saclib_found="no"

fi
done


if test "x$saclib_found" = "xyes" ; then		
	AC_SUBST(SACLIB_CFLAGS)
	AC_SUBST(SACLIB_LIBS)
	AC_DEFINE(HAVE_SACLIB,1,[Define if SACLIB is installed])
	HAVE_SACLIB=yes

	if test "x$saclib_cross" != "xyes"; then
		AC_MSG_RESULT(found)
	else
		AC_MSG_RESULT(unknown)
		echo "WARNING: You appear to be cross compiling, so there is no way to determine"
		echo "whether your SACLIB version is new enough. I am assuming it is."
	fi
	ifelse([$2], , :, [$2])
elif test -n "$saclib_problem"; then
	AC_MSG_RESULT(problem)
	echo "Sorry, your SACLIB version is too old. Disabling."
	ifelse([$3], , :, [$3])
elif test "x$saclib_found" = "xno" ; then
	AC_MSG_RESULT(not found)
	ifelse([$3], , :, [$3])
fi	



AM_CONDITIONAL(LINBOX_HAVE_SACLIB, test "x$HAVE_SACLIB" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])