# Check for LIDIA
# Pascal Giorgi, 2001-12-10
# Modified by Pascal Giorgi, 2003-12-03
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_LIDIA ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for LIDIA Library and define LIDIA_CFLAGS and LIDIA_LIBS

AC_DEFUN([LB_CHECK_LIDIA],
[

AC_ARG_WITH(lidia,
[  --with-lidia=<path>|yes|no Use Lidia library. If argument is no, you do not 
                             have the library installed on your machine (set 
			     as default). If argument is yes or <empty> that 
			     means the library is reachable with the standard
			     search path (/usr or /usr/local). Otherwise you 
			     give the <path> to the directory which contain the 
			     library. 
	     ],
	     [if test "$withval" = yes ; then
			LIDIA_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	      elif test "$withval" != no ; then
			LIDIA_HOME_PATH="$withval"
	     fi],
	     [])

min_lidia_version=ifelse([$1], ,2.1,$1)


dnl Check for existence

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

if test -n "$LIDIA_HOME_PATH" ; then
AC_MSG_CHECKING(for LIDIA >= $min_lidia_version)
fi

for LIDIA_HOME in ${LIDIA_HOME_PATH} 
 do	
if test -r "$LIDIA_HOME/include/LiDIA/LiDIA.h"; then
	if test "x$LIDIA_HOME" != "x/usr" -a "x$LIDIA_HOME" != "x/usr/local"; then
		LIDIA_CFLAGS="-I${LIDIA_HOME}/include"
		LIDIA_LIBS="-L${LIDIA_HOME}/lib -lLiDIA"
	else
		LIDIA_CFLAGS=
		LIDIA_LIBS="-lLiDIA"		
	fi	
	CXXFLAGS="${BACKUP_CXXFLAGS} ${LIDIA_CFLAGS} ${GMP_CFLAGS}" 
	LIBS="${BACKUP_LIBS} ${LIDIA_LIBS} ${GMP_LIBS}"

	AC_TRY_LINK(
	[#include <LiDIA/bigint.h>],
	[LiDIA::bigint a;],
	[
	AC_TRY_RUN(
	[#include <LiDIA/LiDIA.h>
	int main () {  if (LIDIA_MAJOR_VERSION < 2) return -1; else return 0; }
	],[
	lidia_found="yes"
	break
	],[	
	lidia_problem="$problem $LIDIA_HOME"	
	unset LIDIA_CFLAGS
	unset LIDIA_LIBS
	],[
	lidia_found="yes"
	lidia_cross="yes"
	break
	])	
	],
	[
	lidia_found="no"
	lidia_checked="$checked $LIDIA_HOME"
	unset LIDIA_CFLAGS
	unset LIDIA_LIBS	
	])
else
	lidia_found="no"
fi
done

if test "x$lidia_found" = "xyes" ; then		
	AC_SUBST(LIDIA_CFLAGS)
	AC_SUBST(LIDIA_LIBS)
	AC_DEFINE(HAVE_LIDIA,1,[Define if LIDIA is installed])
	HAVE_LIDIA=yes
	if test "x$lidia_cross" != "xyes"; then
		AC_MSG_RESULT(found)
	else
		AC_MSG_RESULT(unknown)
		echo "WARNING: You appear to be cross compiling, so there is no way to determine"
		echo "whether your LIDIA version is new enough. I am assuming it is."
	fi
	ifelse([$2], , :, [$2])	
elif test -n "$lidia_problem"; then
	AC_MSG_RESULT(problem)
	echo "Sorry, your LIDIA version is too old. Disabling."
	ifelse([$3], , :, [$3])
elif test "x$lidia_found" = "xno" ; then	
	AC_MSG_RESULT(not found)
	ifelse([$3], , :, [$3])
fi	


AM_CONDITIONAL(LINBOX_HAVE_LIDIA, test "x$HAVE_LIDIA" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
