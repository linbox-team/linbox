export classes_destination = $(CURDIR)/WEB-INF/classes
special_make = ComputationServer.war
to_clean = $(special_make)

ComputationServer.war : all
	$(JAVADIR)jar cvf ComputationServer.war .

.PHONY : all
all :
	cd JAVA-SRC && $(MAKE)

.PHONY : clean
clean :
	-rm $(to_clean)
	-rm -r $(classes_destination)/*
