# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence

ACLOCAL_AMFLAGS = -I macros

SUBDIRS=linbox macros tests interfaces doc examples
#ACLOCAL_AM_FLAGS=-I m4

docs:doc/linbox-html/index.html

doc/linbox-html/index.html:
	(cd doc; make docs)

examples:
	(cd examples; make examples)

.PHONY:examples

bin_SCRIPTS=linbox-config

VERSION=1.1.7
