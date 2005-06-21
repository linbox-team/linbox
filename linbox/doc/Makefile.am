if LINBOX_BUILD_DOC
USE_TARGETS = doc
INSTALL_TARGETS = install-doc
endif
docdir=$(LINBOX_DOC_PATH)

all all-local: $(USE_TARGETS)

install install-data-local: $(INSTALL_TARGETS)

doc:
	doxygen 
	cp tutorial.html html
	cp install-dist.html html
	cp install-dev.html html
	cp ../INSTALL html
	cp ../COPYING html
	cp ../AUTHORS html
	cp ../HACKING html

install-doc :
	cp -r html $(docdir)/html
	cp linbox.html $(docdir)/linbox.html

EXTRA_DIST=\
	linbox.html		\
	tutorial.html		\
	install-dev.html	\
	install-dist.html	\
	Doxyfile		\
	doc.doxy		\
	mainpage.doxy


clean-local:
	rm -rf html


