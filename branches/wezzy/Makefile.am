# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence

ACLOCAL_AMFLAGS = -I macros

SUBDIRS=linbox macros tests interfaces doc examples benchmarks

# I intend to pare down the target "check" so as to be a smaller check used
# at installation time with a focus on confirming that dependencies on
# other packages and compilation environment are in order.  In contrast, the
# target "fullcheck" should become as exhaustive as possible to serve
# developers for regression testing. -bds 2011Feb ...also benchmarks will be.

#  Fullcheck provides a full check for regression testing purposes.
fullcheck:
	(make ; cd tests; make fullcheck)

docs:doc/linbox-html/index.html

doc/linbox-html/index.html:
	(cd doc; make docs)

docs_dev:doc/linbox-dev-html/pages.html

doc/linbox-dev-html/pages.html:
	(cd doc; make docs_dev)

examples:
	(cd examples; make examples)

benchmarks:
	(cd benchmarks; make benchmarks)

.PHONY:examples benchmarks

bin_SCRIPTS=linbox-config

bundled:
	sh auto-install.sh

uninstall-hook:
	(test -d "$(includedir)/linbox" && rmdir "$(bindir)" \
		"$(libdir)" \
		"$(mandir)/man1" \
		"$(mandir)" \
		"$(includedir)/linbox/algorithms" \
		"$(includedir)/linbox/blackbox" \
		"$(includedir)/linbox/element" \
		"$(includedir)/linbox/field/Givaro" \
		"$(includedir)/linbox/field/Modular" \
		"$(includedir)/linbox/field/NTL" \
		"$(includedir)/linbox/field/" \
		"$(includedir)/linbox/kaapi" \
		"$(includedir)/linbox/matrix" \
		"$(includedir)/linbox/randiter" \
		"$(includedir)/linbox/ring" \
		"$(includedir)/linbox/solutions" \
		"$(includedir)/linbox/switch" \
		"$(includedir)/linbox/util/formats" \
		"$(includedir)/linbox/util" \
		"$(includedir)/linbox/vector" \
		"$(includedir)/linbox" ) || true

VERSION=1.3.0

EXTRA_DIST=auto-install.sh
#incremente-versions
