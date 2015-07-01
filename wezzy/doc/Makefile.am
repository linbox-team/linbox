# Copyright (c) 2010 the LinBox group
#
# ========LICENCE========
# This file is part of the library LinBox.
#
# LinBox is free software: you can redistribute it and/or modify
# it under the terms of the  GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========
#/




if LINBOX_BUILD_DOC
USE_TARGETS = docs
INSTALL_TARGETS = install-doc
endif
docdir=$(LINBOX_DOC_PATH)

man1_MANS = linbox-config.1

all all-local: $(USE_TARGETS)

install install-data-local: $(INSTALL_TARGETS)

docs :
	sed -i 's/\\version\s1.*/\\version\ $(VERSION)/' mainpage.doxy
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
	mkdir -p $(DESTDIR)/$(docdir)
	cp -rp linbox-html $(DESTDIR)/$(docdir)/linbox-html
	cp -rp linbox-dev-html $(DESTDIR)/$(docdir)/linbox-dev-html
	cp -p linbox.html $(DESTDIR)/$(docdir)/linbox.html

EXTRA_DIST=             \
    linbox.html         \
    tutorial.html       \
    install-dev.html    \
    index-dev.html      \
    install-dist.html   \
    Doxyfile.mod        \
    DoxyfileDev.mod     \
    doc.doxy            \
    mainpage.doxy       \
    tutorial.doxy       \
    linbox-config.1


clean-local :
	rm -rf linbox-html
	rm -rf linbox-dev-html


