SUBDIRS=library tools wrappers

EXTRA_DIST=linbox.dxx

# make subdir flat/ which contains a link to every code file
flat: $top_srcdir/macros/flatten.sh
	$top_srcdir/macros/flatten.sh
