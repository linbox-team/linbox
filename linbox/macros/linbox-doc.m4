
AC_DEFUN([LB_DOC],
[

AC_MSG_CHECKING(whether to build documentation)


AC_ARG_WITH(docdir,
	    [ --with-docdir=<path> 
					   Where the LinBox documentation should be installed
	    ],
            [
		LINBOX_DOC_PATH="$withval"
	    ],
	    [
		eval LINBOX_DOC_PATH="${prefix}/doc"
	    ])

AC_SUBST(LINBOX_DOC_PATH)

AC_ARG_WITH(doxygen,
	    [ --with-doxygen=<path> 
					   Give the path to Doxygen. 
					   Note: --enable-doc needed					 
	    ],
            [
		DOXYGEN_PATH="$PATH $withval"
	    ],
	    [
		DOXYGEN_PATH="$PATH"
	    ])

AC_ARG_ENABLE(doc,[--enable-doc Enable building documentation],
[
AC_MSG_RESULT(yes)
AC_CHECK_PROG(USE_DOXYGEN,doxygen,"found","not found", $DOXYGEN_PATH)
if test "x$USE_DOXYGEN" = "xnot found"; then
   echo '*******************************************************************************'
   echo 'Documentation will not be built.'
   echo '*******************************************************************************'
   AM_CONDITIONAL(LINBOX_BUILD_DOC, false)
else
   AM_CONDITIONAL(LINBOX_BUILD_DOC, true)	
fi
],
[
AC_MSG_RESULT(no)
AM_CONDITIONAL(LINBOX_BUILD_DOC, false)
])
])
