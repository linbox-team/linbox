special_make = Locations.java
to_clean = $(special_make) *.class

TIMEOUT = 1

ifndef COMPUTING_SERVER_TEMP_DIR
ifdef TEMP_DIR
COMPUTING_SERVER_TEMP_DIR = $(TEMP_DIR)
endif
endif

tell_to_define = @echo $(1) not defined.  Use \'$(MAKE) $(1)=[value]\' \
                 && exit 1
sed_replace = 's/\$${$(1)}/$(subst /,\/,$(2))/g'

.PHONY : all
all : $(special_make)
	$(JAVADIR)javac *.java

Locations.java : Locations.java.unparsed
ifdef noclobber
	unset noclobber
endif
ifndef PORT
	$(call tell_to_define,PORT)
else
ifndef TIMEOUT
	$(call tell_to_define,TIMEOUT)
else
ifndef LINBOX_FOLDER
	$(call tell_to_define,LINBOX_FOLDER)
else
ifndef COMPUTING_SERVER_TEMP_DIR
	$(call tell_to_define,COMPUTING_SERVER_TEMP_DIR)
else
	sed -e $(call sed_replace,PORT,$(PORT)) \
	    -e $(call sed_replace,TIMEOUT,$(TIMEOUT)) \
	    -e $(call sed_replace,LINBOX_FOLDER,$(LINBOX_FOLDER)) \
	    -e $(call sed_replace,COMPUTING_SERVER_TEMP_DIR,$(COMPUTING_SERVER_TEMP_DIR)) \
	    $< > $@
endif
endif
endif
endif

.PHONY : clean
clean :
	-rm $(to_clean)

%.class : %.java
	$(JAVADIR)javac $<
