# Copyright(c)'2000 The Linbox Group
# Jean-Guillaume.Dumas@imag.fr

EXEC=birthdaytest muple_test

all lib: ${EXEC}

tests:all
	@echo "******* Testing random generators ..."
	@echo ""
	./birthdaytest
	./muple_test
	@echo "****************************************"

header:
clean distclean:
	- \rm -f ${EXEC}
 
include ${LINBOXROOT}/include/linbox_makefile

SPARCFLAGS=-mcpu=ultrasparc -mtune=ultrasparc
OPTFLAGS=-O5 -funroll-all-loops -felide-constructors
LINBOX_OPTFLAGS=$(OPTFLAGS) $(SPARCFLAGS)

