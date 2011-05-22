# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence

ACLOCAL_AMFLAGS = -I macros

SUBDIRS=linbox macros tests interfaces doc examples

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

docs_dev:doc/linbox-dev-html/index.html

doc/linbox-dev-html/index.html:
	(cd doc; make docs_dev)

examples:
	(cd examples; make examples)


.PHONY:examples

bin_SCRIPTS=linbox-config

bundled:
	sh auto-install.sh

VERSION=1.2.0

EXTRA_DIST=auto-install.sh
#incremente-versions
