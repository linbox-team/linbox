# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
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
