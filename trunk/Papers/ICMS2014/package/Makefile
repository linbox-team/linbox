FOLDER=Boyer_Brice
all:
	pdflatex root_$(FOLDER) ;
	pdflatex root_$(FOLDER) ;
	-mv -f root_$(FOLDER).pdf $(FOLDER)/ ;
	-cp -f root_$(FOLDER).tex $(FOLDER)/
clean:
	-rm -f *log *aux	
