SUBDIRS=gmp++ linbox macros tests examples interfaces doc

doc:doc/linbox-html/index.html

doc/linbox-html/index.html:
	(cd doc; make docs)

examples:examples/dot-product

examples/dot-product:linbox-config.h
	(cd examples; make examples)

bin_SCRIPTS=linbox-config

include_HEADERS=linbox-config.h HACKING

VERSION=1.0.6
