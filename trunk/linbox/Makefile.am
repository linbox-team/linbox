SUBDIRS=gmp++ linbox macros tests examples interfaces doc

doc:doc/linbox-html/index.html

doc/linbox-html/index.html:
	(cd doc; make docs)

examples:
	(cd examples; make examples)

bin_SCRIPTS=linbox-config

include_HEADERS=linbox-config.h HACKING

VERSION=0.99.1
