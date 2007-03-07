# Check for GMP
# Modified by Pascal Giorgi, 2003-12-03

dnl LB_CHECK_GMP ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for the GNU Multiprecision library and define GMP_CFLAGS and GMP_LIBS

AC_DEFUN([LB_CHECK_GMP],
[

AC_ARG_WITH(gmp,
		[  --with-gmp= <path>|yes
	   				   Use GMP library. This library is mandatory for LinBox compilation. 
					   If argument is yes or <empty> that means the library is reachable with the standard
					   search path "/usr" or "/usr/local"  (set as default).
	 				   Otherwise you give the <path> to the directory which contain the library. 
		],
		[if test "$withval" = yes ; then
			GMP_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	         elif test "$withval" != no ; then
			GMP_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}"
	        fi],
		[GMP_HOME_PATH="${DEFAULT_CHECKING_PATH}"])

min_gmp_version=ifelse([$1], ,3.1.1,$1)

dnl Check for existence
BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for GMP >= $min_gmp_version)

for GMP_HOME in ${GMP_HOME_PATH} 
  do	
	if test -r "$GMP_HOME/include/gmp.h"; then

		if test "x$GMP_HOME" != "x/usr" -a "x$GMP_HOME" != "x/usr/local"; then
			GMP_CFLAGS="-I${GMP_HOME}/include"
			GMP_LIBS="-L${GMP_HOME}/lib -lgmp"	
		else
			GMP_CFLAGS=
			GMP_LIBS="-lgmp"		
		fi
	
		CXXFLAGS="${CXXFLAGS} ${GMP_CFLAGS}"
		LIBS="${LIBS} ${GMP_LIBS}"

		AC_TRY_LINK(
		[#include <gmp.h>],
		[mpz_t a; mpz_init (a);],
		[
        		AC_TRY_RUN(
 			[#include <gmp.h>
			 int main () {  if (__GNU_MP_VERSION < 3) return -1; else return 0; }
		  	],[
				AC_MSG_RESULT(found)
				AC_SUBST(GMP_CFLAGS)
		  		AC_SUBST(GMP_LIBS)
				AC_DEFINE(HAVE_GMP,1,[Define if GMP is installed])
				# See if we are running GMP 4.0
	   			AC_MSG_CHECKING(whether GMP is 4.0 or greater)
		   		AC_TRY_RUN(
		   		[#include <gmp.h>
	    			int main () { if (__GNU_MP_VERSION < 4) return -1; else return 0; }
	   			],[
					gmp_found="yes"
					AC_MSG_RESULT(yes)
					# See if GMP was compiled with --enable-cxx
					AC_MSG_CHECKING(whether GMP was compiled with --enable-cxx)
					AC_TRY_RUN(
					[#include <gmpxx.h>
					int main () { mpz_class a(2),b(3),c(5); if ( a+b == c ) return 0; else return -1; }
					],[
						AC_MSG_RESULT(yes)
						GMP_VERSION=""
						GMP_LIBS="$GMP_LIBS -lgmpxx"
						AC_SUBST(GMP_VERSION)
					],[
						gmp_found="no"	
						AC_MSG_RESULT(no)
						AC_DEFINE(GMP_NO_CXX,1,[Define if GMP has no <gmpxx.h>])
						unset GMP_CFLAGS
						unset GMP_LIBS	
						GMP_VERSION="-DGMP_NO_CXX"
						AC_SUBST(GMP_VERSION)
					],[
						dnl This should never happen
						AC_MSG_RESULT(no)
					])				
				],[
					AC_MSG_RESULT(no)
					AC_DEFINE(GMP_VERSION_3,1,[Define if GMP is version 3.xxx])
					GMP_VERSION="-DGMP_VERSION_3"
					AC_SUBST(GMP_VERSION)
				],[
					dnl This should never happen
					AC_MSG_RESULT(no)
				])
				ifelse([$2], , :, [$2])
				break
			],[			
				gmp_problem="$gmp_problem $GMP_HOME"
				unset GMP_CFLAGS
				unset GMP_LIBS	
			],[
				AC_MSG_RESULT(unknown)
				echo "WARNING: You appear to be cross compiling, so there is no way to determine"
				echo "whether your GMP version is new enough. I am assuming it is."
				AC_SUBST(GMP_CFLAGS)
				AC_SUBST(GMP_LIBS)
				AC_DEFINE(HAVE_GMP,1,[Define if GMP is installed])	
				ifelse([$2], , :, [$2])
				break
			])	
		],[
		gmp_found="no"	
		unset GMP_CFLAGS
		unset GMP_LIBS	
		])

	else
		gmp_found="no"	
	fi
done

if test "x$gmp_found" != "xyes"; then
	if test -n "$gmp_problem"; then
		AC_MSG_RESULT(problem)
		echo "Sorry, your GMP version is too old. Disabling."
	elif test "x$gmp_found" != "xno"; then
		AC_MSG_RESULT(not found)
	fi
	ifelse($3, , :, $3)
fi


CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
