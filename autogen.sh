#!/bin/sh -ex

srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.

# run autoconf, libtoolize, automake, etc.
autoreconf -vif $srcdir

# run configure script
if test x$NOCONFIGURE = x; then
	$srcdir/configure "$@"
fi
