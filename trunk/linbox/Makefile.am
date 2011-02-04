# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence

ACLOCAL_AMFLAGS = -I macros

SUBDIRS=linbox macros tests interfaces doc examples

# kludge to escape autotools tyranny
#Check: check-recursive
#check-recursive:
#	(cd tests; make tests)

tests: tests/checker
	(cd tests; make tests)

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

VERSION=1.1.7
