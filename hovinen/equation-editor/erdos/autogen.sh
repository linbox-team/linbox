#!/bin/sh
# Run this to generate all the initial makefiles, etc.

srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.
export srcdir

PKG_NAME="Erdos"

(test -f $srcdir/configure.in \
  && test -f $srcdir/erdos.h) || {
    echo -n "**Error**: Directory "\`$srcdir\'" does not look like the"
    echo " top-level erdos directory"
    exit 1
}

. $srcdir/macros/autogen.sh
