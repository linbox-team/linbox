EXTRA_DIST=\
	copy16.gif			\
	tailer_doc++.html		\
	linbox.html			\
	install-dev.html

SUBDIRS = html tex

html: ../linbox/linbox.dxx
	doc++ -d html ../linbox/linbox.dxx

tex: ../linbox/linbox.dxx
	doc++ -t -o doc.tex ../linbox/linbox.dxx

