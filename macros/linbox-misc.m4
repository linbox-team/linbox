# linbox miscellaneous functonnnalities


AC_DEFUN([LB_MISC],
[

AC_ARG_WITH(default,
[  --with-default=<path> Add <path> to the default path for external package 
  		        checking. Set as default with /usr and /usr/local.
],
	    [if test "$withval" = yes ; then
			echo "Default path = /usr /usr/local"
			DEFAULT_CHECKING_PATH="/usr /usr/local"
	      else
			echo "Default path = $withval /usr /usr/local"
			DEFAULT_CHECKING_PATH="$withval /usr /usr/local"
	     fi
	     ],
	     [
		echo "Default path = /usr /usr/local"
		DEFAULT_CHECKING_PATH="/usr /usr/local"
             ])


AC_ARG_WITH(all,
[  --with-all= <path>|yes|no Use all external packages. If the argument is no, 
  	      		   you not sure that all libraries are reachable with 
			   the default path. If the argument is yes or <empty>,
			   that means that all libraries are reachable with the
			   default path. Otherwise add <path> to default path 
			   and enable all external packages.
],
	    [if test "$withval" = yes ; then
			check_all="yes"
			echo "Checking all external packages in ${DEFAULT_CHECKING_PATH}"

	      elif test "$withval" != no ; then
			check_all="yes"
			DEFAULT_CHECKING_PATH="$withval ${DEFAULT_CHECKING_PATH}"
			echo "Checking all external packages in ${DEFAULT_CHECKING_PATH}"			
	     fi
	     ],
	     [])
					
if test -n "$check_all"; then

	GMP_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	GIVARO_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	NTL_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	LIDIA_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	SACLIB_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	MAPLE_HOME_PATH="${DEFAULT_CHECKING_PATH} unknown"
#	EXPAT_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	BLAS_HOME_PATH="${DEFAULT_CHECKING_PATH}"
fi


])
