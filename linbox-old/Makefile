all: $(LINBOXROOT)/README header lib

docs : $(LINBOXROOT)/README

$(LINBOXROOT)/README:doc/README.html doc/.sedfile doc/.iso-ent
	sed -f doc/.sedfile doc/README.html | sed -f doc/.iso-ent > $(LINBOXROOT)/README

LAWS=listings header lib clean distclean

$(LAWS):
	$(MAKE) -C src $@

tests:
	$(MAKE) -C src/examples $@


# What is the current lib name ! 
LIBNAME=$(shell basename `pwd`)
# Will tar only pwd/Makefile, pwd/configure, pwd/doc/*, ... 
TOTAR=Makefile configure doc config src 
SUBTARS=$(TOTAR:%=$(LIBNAME)/%)
# tgz file is built in ../

tar:
	@echo "CVS" > /tmp/tar_tmp.exc
	@echo "ARCHIVES" >> /tmp/tar_tmp.exc
	(cd ..; tar -cvf $(LIBNAME).tar --exclude-from=/tmp/tar_tmp.exc $(SUBTARS); gzip -v $(LIBNAME).tar)
	@\rm -rf  /tmp/tar_tmp.exc
