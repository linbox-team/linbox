EXTRA_DIST=\
	copy16.gif			\
	tailer_doc++.html		\
	linbox.html			\
	install-dev.html

SUBDIRS = tex

html: doc.dxx ../linbox.dxx
	doc++ -d html --filenames-path --no-java-graphs ../linbox.dxx
	csh fix-index

tex: ../linbox/linbox.dxx
	doc++ -t -o doc.tex linbox.dxx

