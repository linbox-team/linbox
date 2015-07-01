dnl Copyright (c) the LinBox group
dnl This file is part of LinBox

 dnl ========LICENCE========
 dnl This file is part of the library LinBox.
 dnl
 dnl LinBox is free software: you can redistribute it and/or modify
 dnl it under the terms of the  GNU Lesser General Public
 dnl License as published by the Free Software Foundation; either
 dnl version 2.1 of the License, or (at your option) any later version.
 dnl
 dnl This library is distributed in the hope that it will be useful,
 dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
 dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 dnl Lesser General Public License for more details.
 dnl
 dnl You should have received a copy of the GNU Lesser General Public
 dnl License along with this library; if not, write to the Free Software
 dnl Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 dnl ========LICENCE========
 dnl

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


AC_MSG_CHECKING(whether dot works)
res=yes;
(dot -V) < /dev/null > /dev/null 2>&1 || res=no
AC_MSG_RESULT([$res])
AS_IF([test $res = yes],
[
sed 's/^HAVE_DOT.*/HAVE_DOT = YES/' doc/Doxyfile.mod > doc/Doxyfile
sed 's/^HAVE_DOT.*/HAVE_DOT = YES/' doc/DoxyfileDev.mod > doc/DoxyfileDev
],
[ cp doc/Doxyfile.mod doc/Doxyfile ;
cp doc/DoxyfileDev.mod doc/DoxyfileDev
])


])
