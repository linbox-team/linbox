#!/bin/sh
# Copyright (c) LinBox
#
# ========LICENCE========
# This file is part of the library LinBox.
#
# LinBox is free software: you can redistribute it and/or modify
# it under the terms of the  GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========
#/


# Run this to generate all the initial makefiles, etc.

# Recover command line, with double-quotes
CMDLINE=""
for arg in "$@"
do
    WHO="`echo $arg | cut -d'=' -f1`"
    WHAT="`echo $arg | cut -s -d'=' -f2`"
    if test "x$WHAT" = "x"; then
    	CMDLINE="$CMDLINE $WHO"
    else
	CMDLINE="$CMDLINE $WHO=\"$WHAT\""
    fi
done

echo  "$0 $CMDLINE" > autogen.status
chmod +x autogen.status

# Starts configuring
srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.

PKG_NAME="Linbox Library"

(test -f $srcdir/configure.ac \
	&& test -f $srcdir/linbox/linbox.doxy) || {
	echo -n "**Error**: Directory "\`$srcdir\'" does not look like the"
	echo " top-level "\`$PKG_NAME\'" directory"
	exit 1
}

ORIGDIR=`pwd`
cd $srcdir
PROJECT=linbox
TEST_TYPE=-f

DIE=0

# Defaults
LIBTOOL=libtool
LIBTOOLIZE=libtoolize

(autoconf --version) < /dev/null > /dev/null 2>&1 || {
	echo
	echo "You must have autoconf installed to compile $PROJECT."
	echo "Download the appropriate package for your distribution,"
	echo "or get the source tarball at ftp://ftp.gnu.org/pub/gnu/"
	DIE=1
}

(automake --version) < /dev/null > /dev/null 2>&1 || {
	echo
	echo "You must have automake installed to compile $PROJECT."
	echo "Download the appropriate package for your distribution,"
	echo "or get the source tarball at ftp://ftp.gnu.org/pub/gnu/"
	DIE=1
}


(grep "^AC_PROG_LIBTOOL" configure.ac >/dev/null) && {
  ($LIBTOOL --version) < /dev/null > /dev/null 2>&1 || {
     echo
     echo "**Error**: You must have \`libtool' installed to compile $PROJECT."
     echo "Download the appropriate package for your distribution,"
     echo "or get the source tarball at ftp://ftp.gnu.org/pub/gnu/"
     DIE=1
  }
}

grep "^AM_GNU_GETTEXT" configure.ac >/dev/null && {
  grep "sed.*POTFILES" $srcdir/configure.ac >/dev/null || \
	(gettext --version) < /dev/null > /dev/null 2>&1 || {
	echo
	echo "**Error**: You must have \`gettext' installed to compile $PROJECT."
	echo "Download the appropriate package for your distribution,"
	echo "or get the source tarball at ftp://ftp.gnu.org/pub/gnu/"
	DIE=1
  }
}

if test "$DIE" -eq 1; then
	exit 1
fi


if test -z "$*"; then
	echo "I am going to run ./configure with no arguments - if you wish "
	echo "to pass any to it, please specify them on the $0 command line."
fi

case $CC in
	*xlc | *xlc\ * | *lcc | *lcc\ *) am_opt=--include-deps;;
esac

for coin in `find . -name configure.ac -print`
do 
	dr=`dirname $coin`
	if test -f $dr/NO-AUTO-GEN; then
		echo skipping $dr -- flagged as no auto-gen
	else
		echo processing $dr
		macrodirs=`sed -n -e 's,AM_ACLOCAL_INCLUDE(\(.*\)),\1,gp' < $coin`
		( cd $dr
		aclocalinclude="$ACLOCAL_FLAGS"
		for k in $macrodirs; do
			if test -d $k; then
				aclocalinclude="$aclocalinclude -I $k"
				##else 
				##  echo "**Warning**: No such directory \`$k'.  Ignored."
			fi
		done
		if grep "^AM_GNU_GETTEXT" configure.ac >/dev/null; then
			if grep "sed.*POTFILES" configure.ac >/dev/null; then
				: do nothing -- we still have an old unmodified configure.ac
			else
				echo "Creating $dr/aclocal.m4 ..."
				test -r $dr/aclocal.m4 || touch $dr/aclocal.m4
				echo "Running gettextize...  Ignore non-fatal messages."
				echo "no" | gettextize --force --copy
				echo "Making $dr/aclocal.m4 writable ..."
				test -r $dr/aclocal.m4 && chmod u+w $dr/aclocal.m4
			fi
		fi
		if grep "^AM_GNOME_GETTEXT" configure.ac >/dev/null; then
			echo "Creating $dr/aclocal.m4 ..."
			test -r $dr/aclocal.m4 || touch $dr/aclocal.m4
			echo "Running gettextize...  Ignore non-fatal messages."
			echo "no" | gettextize --force --copy
			echo "Making $dr/aclocal.m4 writable ..."
			test -r $dr/aclocal.m4 && chmod u+w $dr/aclocal.m4
		fi
		if grep "^AC_PROG_LIBTOOL" configure.ac >/dev/null; then
			echo "Running libtoolize..."
			$LIBTOOLIZE --force --copy
		fi
		echo "Running aclocal $aclocalinclude ..."
		aclocal $aclocalinclude
		if grep "^AC_CONFIG_HEADERS" configure.ac >/dev/null; then
			echo "Running autoheader..."
			autoheader
		fi
		echo "Running automake --gnu $am_opt ..."
		automake -c --add-missing --gnu $am_opt
		echo "Running autoconf ..."
		autoconf
		)
	fi
done

conf_flags="--enable-maintainer-mode" 
#--enable-iso-c

cd "$ORIGDIR"

if test x$NOCONFIGURE = x; then
	echo Running $srcdir/configure $conf_flags "$@" ...
	$srcdir/configure $conf_flags "$@" \
		&& echo "Now type \`make' to compile $PROJECT"  || exit 1
else
	echo Skipping configure process.
fi

