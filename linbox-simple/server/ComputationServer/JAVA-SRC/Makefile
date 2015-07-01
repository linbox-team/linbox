SUBDIRS = LinBoxServerInterface
special_make =
to_clean = $(special_make) *.class

.PHONY : all
all :
	$(foreach dir,$(SUBDIRS),cd $(dir) && $(MAKE) &&) \
	$(JAVADIR)javac -d $(classes_destination) *.java

.PHONY : clean
clean :
	-rm $(to_clean)

%.class : %.java
	$(JAVADIR)javac $<
