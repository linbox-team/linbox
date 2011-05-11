# Copyright (c) 2011 the LinBox group
# Brice Boyer <bboyer@imag.fr>
# Adapted from linbox-doc.m4
# This file is part of LinBox
# see COPYING for licence


# Copyright (c) the LinBox group
# This file is part of LinBox
# see COPYING for licence

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

