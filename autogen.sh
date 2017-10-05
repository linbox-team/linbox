#!/bin/sh -ex

srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.

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

echo  "$0 $CMDLINE" > $srcdir/autogen.status
chmod +x $srcdir/autogen.status

# run autoconf, libtoolize, automake, etc.
autoreconf -vif $srcdir

# run configure script
if test x$NOCONFIGURE = x; then
	$srcdir/configure "$@"
fi
