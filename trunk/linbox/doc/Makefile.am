EXTRA_DIST=\
	copy16.gif			\
	tailer_doc++.html		\
	linbox.html			\
	install-dev.html

SUBDIRS = tex

doc: 
	rm -r html
	gmake html

html: doc.dxx ../linbox.dxx
	doc++ -w -d html --filenames-path --no-java-graphs ../linbox.dxx
	csh fix-index

doc.tex: ../linbox.dxx
	doc++ -t -o doc.tex ../linbox.dxx

