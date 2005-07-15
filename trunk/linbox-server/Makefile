.PHONY : all
all : computing-server, server-interface

.PHONY : computing-server
computing-server :
	cd LinBoxComputingServer && $(MAKE)

.PHONY : server-interface
server-interface :
	cd ComputationServer && $(MAKE)

%.class : %.java
	javac $<

.DELETE_ON_ERROR :
