if USE_DOCPLUSPLUS_TEX
USE_TARGETS = linbox.tex linbox.ps linbox.dvi
endif

all-local: $(USE_TARGETS)

linbox.ps: linbox.dvi
	dvips linbox.dvi -o linbox.ps

linbox.dvi: linbox.tex
	latex "\scrollmode\input linbox"; exit 0

linbox.tex: $(top_srcdir)/linbox/linbox.dxx
	doc++ -t -o linbox.tex $(top_srcdir)/linbox/linbox.dxx

