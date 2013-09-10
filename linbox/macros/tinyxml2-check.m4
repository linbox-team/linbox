dnl Check for FPLLL
dnl Copyright (c) the LinBox group
dnl Boyer Brice <bbboyer@ncsu.edu> 10/09/13
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

AC_DEFUN([LB_CHECK_XML],
[
AC_MSG_CHECKING(if tinyxml2 is available)
SAVED_LIBS=$LIBS
    LIBS="$LIBS -ltinyxml2"
    AC_TRY_LINK(
      [#include <tinyxml2.h>],
      [
      using namespace tinyxml2 ;
      XMLDocument doc;
      ],
      [
AC_MSG_RESULT( OK.)
AC_DEFINE(HAVE_TINYXML2,1,[Define if tinyxml2 is installed])
XML_LIBS="-ltinyxml2"
AC_SUBST(XML_LIBS)
],
      [AC_MSG_WARN([tinyxml2 is not installed (no import/export of benchmarks).])]
    )
    LIBS=$SAVED_LIBS
])
