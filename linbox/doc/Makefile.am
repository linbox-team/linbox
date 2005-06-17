if LINBOX_BUILD_DOC

USE_TARGETS = doc
docdir=$(LINBOX_DOC_PATH)

all-local: $(USE_TARGETS)

doc:
	doxygen 
	cp tutorial.html html
	cp install-dist.html html
	cp install-dev.html html
	cp ../INSTALL html
	cp ../COPYING html
	cp ../AUTHORS html
	cp ../HACKING html

EXTRA_DIST=\
	linbox.html		\
	install-dev.html	\
	install-dist.html	\
	Doxyfile		\
	doc.doxy		\
	mainpage.doxy


install-data-local:
	cp -r html $(docdir)/html
	cp linbox.html $(docdir)/linbox.html

clean-local:
	rm -rf html


endif
