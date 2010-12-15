if LINBOX_BUILD_DOC
USE_TARGETS = docs
INSTALL_TARGETS = install-doc
endif
docdir=$(LINBOX_DOC_PATH)

man1_MANS = linbox-config.1

all all-local: $(USE_TARGETS)

install install-data-local: $(INSTALL_TARGETS)

docs :
	sed -i 's/^\\version.*/\\version\ $(VERSION)/' mainpage.doxy
	if test -d linbox-html ; then echo exists; else mkdir linbox-html ; fi
	if test -d linbox-dev-html ; then echo exists; else mkdir linbox-dev-html ; fi
	cp index-dev.html linbox-dev-html/index.html
	cp tutorial.html linbox-html/
	cp install-dist.html linbox-html/
	cp install-dev.html linbox-html/
	cp ../INSTALL linbox-html/
	cp ../COPYING linbox-html/
	cp ../AUTHORS linbox-html/
	doxygen Doxyfile

docs_dev :
	make docs
	cp tutorial.html linbox-dev-html/
	cp install-dist.html linbox-dev-html/
	cp install-dev.html linbox-dev-html/
	cp ../INSTALL linbox-dev-html/
	cp ../COPYING linbox-dev-html/
	cp ../AUTHORS linbox-dev-html/
	doxygen DoxyfileDev



install-doc :
	mkdir $(docdir)
	cp -r linbox-html $(docdir)/linbox-html
	cp -r linbox-dev-html $(docdir)/linbox-dev-html
	cp linbox.html $(docdir)/linbox.html

EXTRA_DIST=             \
    linbox.html         \
    tutorial.html       \
    install-dev.html    \
    index-dev.html      \
    install-dist.html   \
    Doxyfile            \
    DoxyfileDev         \
    doc.doxy            \
    mainpage.doxy       \
    tutorial.doxy       \
    linbox-config.1


clean-local :
	rm -rf linbox-html
	rm -rf linbox-dev-html


