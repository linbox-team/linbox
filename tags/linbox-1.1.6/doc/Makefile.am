if LINBOX_BUILD_DOC
USE_TARGETS = docs
INSTALL_TARGETS = install-doc
endif
docdir=$(LINBOX_DOC_PATH)

man1_MANS = linbox-config.1 

all all-local: $(USE_TARGETS)

install install-data-local: $(INSTALL_TARGETS)

docs: 
	doxygen 
	cp tutorial.html linbox-html
	cp install-dist.html linbox-html
	cp install-dev.html linbox-html
	cp ../INSTALL linbox-html
	cp ../COPYING linbox-html
	cp ../AUTHORS linbox-html

install-doc :
	mkdir $(docdir)
	cp -r linbox-html $(docdir)/linbox-html
	cp linbox.html $(docdir)/linbox.html

EXTRA_DIST=\
	linbox.html		\
	tutorial.html		\
	install-dev.html	\
	install-dist.html	\
	Doxyfile		\
	doc.doxy		\
	mainpage.doxy		\
	linbox-config.1 


clean-local:
	rm -rf linbox-html


