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



ACLOCAL_AMFLAGS = -I macros

SUBDIRS=linbox macros tests interfaces doc examples benchmarks

# I intend to pare down the target "check" so as to be a smaller check used
# at installation time with a focus on confirming that dependencies on
# other packages and compilation environment are in order.  In contrast, the
# target "fullcheck" should become as exhaustive as possible to serve
# developers for regression testing. -bds 2011Feb ...also benchmarks will be.

#  Fullcheck provides a full check for regression testing purposes.
fullcheck:
	(${MAKE} ; cd tests; ${MAKE} fullcheck)

docs:doc/linbox-html/index.html

doc/linbox-html/index.html:
	(cd doc; ${MAKE} docs)

docs_dev:doc/linbox-dev-html/pages.html

doc/linbox-dev-html/pages.html:
	(cd doc; ${MAKE} docs_dev)

examples:
	(  ${MAKE} ; cd examples; ${MAKE} examples)

benchmarks:
	(cd benchmarks; ${MAKE} benchmarks)

.PHONY:examples benchmarks

bin_SCRIPTS=linbox-config

bundled:
	sh auto-install.sh

uninstall-hook:
	(test -d "$(includedir)/linbox" && rmdir "$(bindir)" \
		"$(libdir)" \
		"$(mandir)/man1" \
		"$(mandir)" \
		"$(includedir)/linbox/algorithms/gauss" \
		"$(includedir)/linbox/algorithms" \
		"$(includedir)/linbox/blackbox" \
		"$(includedir)/linbox/element" \
		"$(includedir)/linbox/field/Givaro" \
		"$(includedir)/linbox/field/Modular" \
		"$(includedir)/linbox/field/NTL" \
		"$(includedir)/linbox/field/" \
		"$(includedir)/linbox/kaapi" \
		"$(includedir)/linbox/matrix" \
		"$(includedir)/linbox/randiter" \
		"$(includedir)/linbox/ring" \
		"$(includedir)/linbox/solutions" \
		"$(includedir)/linbox/switch" \
		"$(includedir)/linbox/util/formats" \
		"$(includedir)/linbox/util" \
		"$(includedir)/linbox/vector" \
		"$(includedir)/linbox" ) || true

VERSION=1.3.1

EXTRA_DIST=auto-install.sh
#incremente-versions
