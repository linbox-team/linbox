# $Source$
# Copyright(c)'94-97 by Givaro Team
# Author: T. Gautier
# $Id$
# Description:
# This file is the makefile to build test examples of
# the Givaro library.

EXEC=rank_sparse_reord_gauss

all:${EXEC}

tests:all
	@echo "******* Testing gaussian elimination ..."
	@echo ""
	./rank_sparse_reord_gauss Examples/cis.mk9.b3.945x1260.sms.gz 32749
	./rank_sparse_reord_gauss Examples/cis.mk9.b3.945x1260.sms.gz 3
	@echo "****************************************"

clean distclean:
	- \rm -f ${EXEC}

include ${GIVAROROOT}/include/givaro_makefile
include ${LINBOXROOT}/include/linbox_makefile

OPTFLAGS=-O5 -funroll-all-loops -felide-constructors -mcpu=ultrasparc -mtune=ultrasparc
OPT+=$(OPTFLAGS)
