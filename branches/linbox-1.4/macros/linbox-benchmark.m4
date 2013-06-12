dnl Copyright (c) 2011 the LinBox group
dnl Brice Boyer <bboyer@imag.fr>
dnl Adapted from linbox-doc.m4
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

AC_DEFUN([LB_BENCH],
[

dnl  AC_MSG_CHECKING(whether to try and benchmark linbox\n)

AC_ARG_WITH(benchdir,
[AC_HELP_STRING([--with-benchdir=<path>], [Where the LinBox benchmarks should be performed])],
            [
		LINBOX_BENCH_PATH="$withval"
	    ],
	    [
		eval LINBOX_BENCH_PATH="${prefix}/benchmarks"
	    ])

AC_SUBST(LINBOX_BENCH_PATH)

AC_ARG_WITH(gnuplot,
[AC_HELP_STRING([--with-gnuplot=<path>], [Give the path to Gnuplot. ])],
            [
		GNUPLOT_PATH="$PATH $withval"
	    ],
	    [
		GNUPLOT_PATH="$PATH"
	    ])
AC_ARG_WITH(ghostscript,
[AC_HELP_STRING([--with-ghostscript=<path>], [Give the path to ghostscript. ])],
            [
		GHOSTSCRIPT_PATH="$PATH $withval"
	    ],
	    [
		GHOSTSCRIPT_PATH="$PATH"
	    ])

dnl  AC_ARG_ENABLE(benchmarks,[ --enable-benchmarks Enables benchmarking],
dnl  [
dnl  AC_MSG_RESULT(yes)
AC_MSG_CHECKING(whether gnuplot works)
res=yes;
export PATH=$GNUPLOT_PATH
(gnuplot --version) < /dev/null > /dev/null 2>&1 || res=no
AC_MSG_RESULT([$res])
if test $res = no ; then
	echo
	echo "You must have gnuplot installed to create benchmark  "
	echo "graphics for LinBox. Download the appropriate package"
	echo "for your distribution, or get the source tarball from"
    echo "http://www.gnuplot.info/download.html                "
else
AC_DEFINE(HAVE_GNUPLOT, 1, [gnuplot available as external program])
fi

AC_MSG_CHECKING(whether ps2pdf works)
res=yes;
export PATH=$GHOSTSCRIPT_PATH
(ps2pdf --version -) < /dev/null > /dev/null 2>&1 || res=no
AC_MSG_RESULT([$res])
if test $res = no ; then
	echo
	echo "You must have ps2pdf installed to create pdf benchmarks"
	echo "graphics for LinBox. Download the appropriate package  "
	echo "for your distribution, or get the source tarball from  "
    echo "http://pages.cs.wisc.edu/~ghost/                       "
else
AC_DEFINE(HAVE_GHOSTSCRIPT, 1, [ps2pdf available as external program])
fi

dnl  AM_CONDITIONAL(LINBOX_BUILD_BENCH, true)
dnl  ],
dnl  [
dnl  AC_MSG_RESULT(no)
dnl  AM_CONDITIONAL(LINBOX_BUILD_BENCH, false)
dnl  ])
])

