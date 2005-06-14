if LINBOX_BUILD_DOC

USE_TARGETS = linbox-html
docdir=$(LINBOX_DOC_PATH)

all-local: $(USE_TARGETS)

linbox-html:
	doxygen 

EXTRA_DIST=\
	linbox.html	\
	install-dev.html	\
	Doxyfile	\
	doc.doxy	\
	mainpage.doxy


install-data-local:
	cp -r html $(docdir)/html
	cp linbox.html $(docdir)/linbox.html

clean-local:
	rm -rf html


endif
