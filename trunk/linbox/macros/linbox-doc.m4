dnl Copyright (c) the LinBox group
dnl This file is part of LinBox
dnl see COPYING for licence

AC_DEFUN([LB_DOC],
[

AC_MSG_CHECKING(whether to build documentation)


AC_ARG_WITH(docdir,
[AC_HELP_STRING([--with-docdir=<path>], [Where the LinBox documentation should be installed])],
            [
		LINBOX_DOC_PATH="$withval"
	    ],
	    [
		eval LINBOX_DOC_PATH="${prefix}/doc"
	    ])

AC_SUBST(LINBOX_DOC_PATH)

AC_ARG_WITH(doxygen,
[AC_HELP_STRING([--with-doxygen=<path>], [Give the path to Doxygen. Note: --enable-doc needed])],
            [
		DOXYGEN_PATH="$PATH $withval"
	    ],
	    [
		DOXYGEN_PATH="$PATH"
	    ])

AC_ARG_ENABLE(doc,[AC_HELP_STRING([--enable-doc], [Enable building documentation])])
WANTDOC="no"
AS_IF([ test "x$enable_doc" = "xyes"],[
   AC_MSG_RESULT(yes)
   AC_MSG_CHECKING(whether doxygen works)
   export PATH=$DOXYGEN_PATH
   (doxygen --version) < /dev/null > /dev/null 2>&1 || {
       AC_MSG_RESULT(no)
       echo
       echo "You must have doxygen installed to create documentation for"
       echo "LinBox. This error only happens if you use --enable-doc."
       echo "Download the appropriate package for your distribution, or get"
       echo "the source tarball from http://www.stack.nl/~dimitri/doxygen/"
       exit -1
   }
   AC_MSG_RESULT(yes)
   WANTDOC="yes"
],[AC_MSG_RESULT(no)])

AM_CONDITIONAL(LINBOX_BUILD_DOC, test "x$WANTDOC" != "xno" )
])
