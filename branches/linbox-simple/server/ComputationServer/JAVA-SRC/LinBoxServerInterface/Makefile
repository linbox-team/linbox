special_make = Locations.java
to_clean = $(special_make) *.class

HOST_NAME = localhost
BASE_ADDRESS = http://localhost:8080/ComputationServer

ifndef MATRIX_DIR
ifdef TEMP_DIR
MATRIX_DIR = $(TEMP_DIR)
endif
endif

ifndef COMPUTATION_DIR
ifdef TEMP_DIR
COMPUTATION_DIR = $(TEMP_DIR)
endif
endif

tell_to_define = @echo $(1) not defined.  Use \'$(MAKE) $(1)=[value]\' \
                 && exit 1
sed_replace = 's/\$${$(1)}/$(subst /,\/,$(2))/g'

.PHONY : all
all : $(special_make)
	$(JAVADIR)javac $(JAVAFLAGS) -d $(classes_destination) *.java

Locations.java : Locations.java.unparsed
ifdef noclobber
	unset noclobber
endif
ifndef PORT
	$(call tell_to_define,PORT)
else
ifndef HOST_NAME
	$(call tell_to_define,HOST_NAME)
else
ifndef BASE_ADDRESS
	$(call tell_to_define,BASE_ADDRESS)
else
ifndef MATRIX_DIR
	$(call tell_to_define,MATRIX_DIR)
else
ifndef COMPUTATION_DIR
	$(call tell_to_define,COMPUTATION_DIR)
else
	sed -e $(call sed_replace,PORT,$(PORT)) \
	    -e $(call sed_replace,HOST_NAME,$(HOST_NAME)) \
	    -e $(call sed_replace,BASE_ADDRESS,$(BASE_ADDRESS)) \
	    -e $(call sed_replace,MATRIX_DIR,$(MATRIX_DIR)) \
	    -e $(call sed_replace,COMPUTATION_DIR,$(COMPUTATION_DIR)) \
	    $< > $@
endif
endif
endif
endif
endif

.PHONY : clean
clean:
	-rm $(to_clean)

%.class : %.java
	$(JAVADIR)javac $(JAVAFLAGS) $<
