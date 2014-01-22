#  TEX=pdflatex
TEX=lualatex
BIB=bibtex
ART=14_linbox

RERUN='(There were undefined references|Rerun to get (bibliographical references|cross-references|the bars|outlines) right)'
UNDEFINED='((Reference|Citation).*undefined)|(Label.*multiply defined)'
#  RERUNBIB = "No file.*\.bbl|Citation.*undefined"
#  #  MAKEIDX = "^[^%]*\\makeindex"
ECHO=/bin/echo



all:
	${TEX} ${ART}.tex
	bibtex ${ART}
	@echo `egrep $(UNDEFINED) ${ART}.log | wc -l` > FAIL
	@if egrep -q $(RERUN) ${ART}.log ; then echo "passe LaTeX" \
		&& ${TEX} ${ART}.tex;  else echo "OK" ; fi
	@echo `egrep $(UNDEFINED) ${ART}.log | wc -l` > FAIL2
	@if [ `cat FAIL` -ne `cat FAIL2` ] ; then echo "passe LaTeX" \
		&& ${TEX} ${ART}.tex; else echo "OK" ; fi
	@echo "=================================================="
	@echo "Citations ou références indéfinies:"
	@egrep -i $(UNDEFINED) ${ART}.log || echo "Aucune"
	@echo "=================================================="
	@echo "Warnings:"
	@egrep --color=auto --line-number Warn ${ART}.log || echo "Aucune"
	@echo "=================================================="
	@echo "Overfull"
	@grep --color=auto --line-number Overfull ${ART}.log || echo "Aucune"
	@echo "=================================================="
	@echo "Underfull"
	@grep --color=auto --line-number Underfull ${ART}.log  || echo "Aucune"
	@echo "=================================================="
	@echo "Completed at: $(shell date)"

once:
	${TEX} ${ART}.tex


redo:
	${TEX} ${ART}.tex
	${TEX} ${ART}.tex
	${BIB} ${ART}
	${TEX} ${ART}.tex

clean :
	-rm -f *.but *.aux  *.blg *.brf *.loa *.lof *.log *.lot *.out *.toc *.idx *.ilg *.ind *.mtc[0-9]* *.maf *.dvi *.mtc *.bmt *.mlf[0-9]* *.mlt[0-9]* *.lol FAIL* *.nav *.snm *.vrb

