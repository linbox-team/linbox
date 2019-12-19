# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# ========LICENCE========
# This file is part of the library LinBox.
#
# LinBox is free software: you can redistribute it and/or modify
# it under the terms of the  GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========
#/


#FLAG128=-D__USE_128bits -D__USE_SIMD=0


DEFCPPFLAGS  = @DEFAULT_CFLAGS@   -DDISABLE_COMMENTATOR
AM_CPPFLAGS = -I$(top_srcdir) $(DEFCPPFLAGS) $(DEPS_CFLAGS)
AM_CPPFLAGS += $(OPTFLAGS) -I../include $(FLAG128) -DLINBOX_DEBUG -DGIVARO_DEBUG -DFFLASFFPACK_DEBUG -UNDEBUG -DDEBUG
LDADD = $(DEPS_LIBS) $(LDFLAGS)
LDADD += $(top_builddir)/linbox/liblinbox.la
LDADD += $(OPTLIBES)

PROGS=evaluation_build_kernel minpoly_minpoly verif_g0kernel evaluation_evaluation\
sequence_decomp_kernel verif_par_verif evaluation_finalchecks\
sequence_sequence verif_precompute evaluation_g0kernel sequence_transpose\
verif_verif_minpoly

EXTRA_PROGRAMS= $(PROGS) 

CLEANFILES=$(EXTRA_PROGRAMS)


programs:  ${EXTRA_PROGRAMS}

evaluation_build_kernel_SOURCES=  ../src/evaluation_build_kernel.cpp 
minpoly_minpoly_SOURCES= 	    ../src/minpoly_minpoly.cpp	    
verif_g0kernel_SOURCES= 	    ../src/verif_g0kernel.cpp	    
evaluation_evaluation_SOURCES=    ../src/evaluation_evaluation.cpp   
sequence_decomp_kernel_SOURCES=   ../src/sequence_decomp_kernel.cpp  
verif_par_verif_SOURCES= 	    ../src/verif_par_verif.cpp	    
evaluation_finalchecks_SOURCES=   ../src/evaluation_finalchecks.cpp  
sequence_sequence_SOURCES= 	    ../src/sequence_sequence.cpp	    
verif_precompute_SOURCES= 	    ../src/verif_precompute.cpp	    
evaluation_g0kernel_SOURCES= 	    ../src/evaluation_g0kernel.cpp	    
sequence_transpose_SOURCES= 	    ../src/sequence_transpose.cpp	    
verif_verif_minpoly_SOURCES=      ../src/verif_verif_minpoly.cpp     

LINBOX=@prefix@

LINBOX_BIN=@bindir@

