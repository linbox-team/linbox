# Check for ATLAS
# Pascal Giorgi , 2003-03-04
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_ATLAS ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for ATLAS and define ATLAS_CFLAGS and ATLAS_LIBS

AC_DEFUN([LB_CHECK_ATLAS],
[

AC_ARG_WITH(atlas,
	    [  --with-atlas=<path>|yes|no 
					   Use Atlas library. 
					   If argument is no, you do not have the library installed on your machine (set as default).
					   If argument is yes or <empty> that means the library is reachable with the standard
					   search path (/usr or /usr/local).
	 				   Otherwise you give the <path> to the directory which contain the library. 	       		
	     ],
	     [if test "$withval" = yes ; then
			ATLAS_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	      elif test "$withval" != no ; then
			ATLAS_HOME_PATH="$withval"
	     fi],
	     [])

AC_ARG_WITH(atlas-includes,
	    [  --with-atlas-includes=<path>
					   Useful for multi-platform installation.
					   Add the path of ATLAS platform dependant headers.					   
					   Need --with-atlas=<generic path> option.
	     ],
	     [ if test -n "$withval" ; then
			OPT_ATLAS_CFLAGS="$withval"
	       fi
	     ],
	     [])

AC_ARG_WITH(atlas-libraries,
	    [  --with-atlas-libraries=<path>
					   Useful for multi-platform installation.
					   Add the path of ATLAS platform dependant libraries.					   
					   Need --with-atlas=<generic path> option.
	     ],
	     [ if test -n "$withval" ; then
			OPT_ATLAS_LIBS="$withval"
	       fi
	     ],
	     [])



min_atlas_version=ifelse([$1], ,3.0,$1)

dnl Check for existence

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

if test -n "$ATLAS_HOME_PATH" ; then
AC_MSG_CHECKING(for ATLAS >= $min_atlas_version)
fi

for ATLAS_HOME in ${ATLAS_HOME_PATH} 
 do
if test -r "$ATLAS_HOME/include/cblas.h"; then

	if test "x$ATLAS_HOME" != "x/usr" -a "x$ATLAS_HOME" != "x/usr/local"; then
		if test -z "${OPT_ATLAS_CFLAGS}" ; then
			ATLAS_CFLAGS="-I${ATLAS_HOME}/include "
		else
			ATLAS_CFLAGS="-I${ATLAS_HOME}/include -I${OPT_ATLAS_CFLAGS}"
		fi
		if test -z "${OPT_ATLAS_LIBS}" ; then
			ATLAS_LIBS="-L${ATLAS_HOME}/lib -llapack -lcblas -latlas" 
		else	
			ATLAS_LIBS="-L${ATLAS_HOME}/lib -L${OPT_ATLAS_LIBS} -llapack -lcblas -latlas"
		fi
		
	else
		if test -z "${OPT_ATLAS_CFLAGS}" ; then
			ATLAS_CFLAGS=
		else
			ATLAS_CFLAGS="-I${OPT_ATLAS_CFLAGS}"
		fi
		if test -z "${OPT_ATLAS_LIBS}" ; then
			ATLAS_LIBS="-llapack -lcblas -latlas" 
		else	
			ATLAS_LIBS="-L${OPT_ATLAS_LIBS} -llapack -lcblas -latlas"
		fi
	fi

	CXXFLAGS="${BACKUP_CXXFLAGS} ${ATLAS_CFLAGS} ${GMP_CFLAGS}" 
	LIBS="${BACKUP_LIBS} ${ATLAS_LIBS} ${GMP_LIBS}"

	AC_TRY_LINK(
	[#include <cblas.h>],
	[double a;],
	[
	AC_TRY_RUN(
	[#include <atlas_buildinfo.h>
	 int main () {  if  (ATL_VERS[0] <3) return -1; else return 0; }
	],[	
	atlas_found="yes"	
	break
	],[	
	atlas_problem="$problem $ATLAS_HOME"	
	unset ATLAS_CFLAGS
	unset ATLAS_LIBS	
	],[
	atlas_found="yes"
	atlas_cross="yes"
	break
	])	
	],
	[
	atlas_found="no"
	atlas_checked="$checked $ATLAS_HOME"
	unset ATLAS_CFLAGS
	unset ATLAS_LIBS
	ifelse([$3], , :, [$3])
	])
else
	atlas_found="no"
fi
done

if test "x$atlas_found" = "xyes" ; then		
	AC_SUBST(ATLAS_CFLAGS)
	AC_SUBST(ATLAS_LIBS)
	AC_DEFINE(HAVE_ATLAS,1,[Define if ATLAS is installed])
	AC_DEFINE(BLAS_AVAILABLE,,[Define if BLAS routines are available])
	HAVE_ATLAS=yes
	if test "x$atlas_cross" != "xyes"; then
		AC_MSG_RESULT(found)
	else
		AC_MSG_RESULT(unknown)
		echo "WARNING: You appear to be cross compiling, so there is no way to determine"
		echo "whether your ATLAS version is new enough. I am assuming it is."
	fi
	ifelse([$2], , :, [$2])
elif test -n "$atlas_problem"; then
	AC_MSG_RESULT(problem)
	echo "Sorry, your ATLAS version is too old. Disabling."
	ifelse([$3], , :, [$3])
elif test "x$atlas_found" = "xno" ; then	
	AC_MSG_RESULT(not found)
	ifelse([$3], , :, [$3])
fi	

AM_CONDITIONAL(LINBOX_HAVE_ATLAS, test "x$HAVE_ATLAS" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])


