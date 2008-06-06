SUBDIRS=gmp++ linbox macros tests interfaces doc examples

docs:doc/linbox-html/index.html

doc/linbox-html/index.html:
	(cd doc; make docs)

examples:examples/dot-product

#examples/dot-product:linbox/linbox-config.h
#	(cd examples; make examples)

bin_SCRIPTS=linbox-config

VERSION=1.1.6rc0
