SUBDIRS = data

AM_CPPFLAGS=-I$(top_srcdir) $(GMP_CFLAGS) $(NTL_CFLAGS) 
LDADD = $(GMP_LIBS) $(NTL_LIBS)

TESTS = test
noinst_PROGRAMS = test

test_SOURCES = test.C
#	test-rational-solver-adaptive.C \
#	test-common.h

#test_smith_form_adaptive_SOURCES=    \
#	test-common.h	\
#	test-smith-form-adaptive.C
