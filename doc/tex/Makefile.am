# Copyright (c) 2010 the LinBox group
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



if USE_DOCPLUSPLUS_TEX
USE_TARGETS = linbox.tex linbox.ps linbox.dvi
endif

CLEANFILES = linbox.tex linbox.ps
MOSTLYCLEANFILES = linbox.aux linbox.log linbox.dvi

all-local: $(USE_TARGETS)

linbox.ps: linbox.dvi
	dvips linbox.dvi -o linbox.ps

linbox.dvi: linbox.tex
	latex "\scrollmode\input linbox"; exit 0

linbox.tex: $(top_srcdir)/linbox/linbox.dxx
	doc++ -t -o linbox.tex $(top_srcdir)/linbox/linbox.dxx
