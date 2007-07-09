SUBDIRS=gmp++ linbox macros tests interfaces doc

docs:doc/linbox-html/index.html

doc/linbox-html/index.html:
	(cd doc; make docs)

examples:examples/dot-product

#examples/dot-product:linbox-config.h
#	(cd examples; make examples)

bin_SCRIPTS=linbox-config

include_HEADERS= HACKING

VERSION=1.1.3
