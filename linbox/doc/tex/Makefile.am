all-local: linbox.tex linbox.ps linbox.dvi

linbox.ps: linbox.dvi
	dvips linbox.dvi -o linbox.ps

linbox.dvi: linbox.tex
	latex -interaction=nonstopmode linbox

linbox.tex: $(top_srcdir)/linbox/linbox.dxx
	doc++ -t -o linbox.tex $(top_srcdir)/linbox/linbox.dxx

