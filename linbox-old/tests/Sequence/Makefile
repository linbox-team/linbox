# $Source$
# Copyright(c)'94-97 by Givaro Team
# Author: T. Gautier
# $Id$
# Description:
# This file is the makefile to build test examples of
# the Givaro library.

EXEC=mas_gfdom mas_zpdom

all:${EXEC}

tests:all 
	@echo "******* Testing Wiedemann's algorithm with two runs and a ppcm ..."
	@echo ""
	./mas_gfdom Examples/cis.mk9.b3.945x1260.sms.gz 3 9 20
	./mas_zpdom Examples/cis.mk9.b3.945x1260.sms.gz 32749 1 20
	@echo "******************************************************************"

clean distclean:
	- \rm -f ${EXEC}
 
include ${GIVAROROOT}/include/givaro_makefile

OPTFLAGS=-O5 -funroll-all-loops -felide-constructors -mcpu=ultrasparc -mtune=ultrasparc -ffast-math -fstrict-aliasing -frerun-loop-opt -fexpensive-optimizations
CXXFLAGS+=-I${HOME}/linbox-1.1/include
OPT+=$(OPTFLAGS)
