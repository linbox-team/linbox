SUBDIRS=library tools wrappers

EXTRA_DIST=linbox.dxx

# make subdir flat/ which contains a link to every code file
flat: $srcdir/macros/flatten.sh
	flatten.sh
