BASE = psd
SRC = ${BASE}.tex

up: psd.ps  
dvi : psd.dvi
ggv: psd.ps ; ggv psd.ps &

${BASE}: ${BASE}.ps
%.dvi:%.tex  *.tex
	latex $*; bibtex $*; latex $*; latex $* 
%.ps:%.dvi
	dvips -P cmz -o $*.ps $*.dvi

${BASE}.dvi:  ${SRC}
	latex ${BASE}; bibtex ${BASE}; latex ${BASE}; latex ${BASE}; 

clean:
	/bin/rm -f *.aux *.log *.blg *.dvi ${BASE}.ps ${BASE}.pdf 

FIGS = 
BIBL = 
STYL = 


TARD = 
tar: 
	- @ mkdir /tmp/$(TARD)
	cp -rf Makefile $(SRC) $(STYL) $(BIBL) $(FIGS) /tmp/$(TARD)/
	@ (cd /tmp/; tar -cvf $(TARD).tar $(TARD))
	@ \rm -rf /tmp/$(TARD)
	@ (cd /tmp/; gzip --best -v $(TARD).tar)
	@ \mv -i /tmp/$(TARD).tar.gz $(TARD).tgz

ispell :
	ispell psd.tex
	ispell abstract.tex
	ispell liegroups.tex
	ispell cra.tex
	ispell application.tex
