INCLUDES=-I$(top_srcdir)

SUBDIRS = data

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) 
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 

TESTS =				\
	test-unparametric-field	\
	test-ntl-RR		\
	test-ntl-ZZ_p		\
	test-ntl-zz_p		\
	test-large-modular	\
	test-vector-domain	\
	test-dense-matrix	\
	test-sparse0	\
	test-diagonal		\
	test-minpoly		\
	test-det		\
	test-rank		\
	test-inverse

noinst_PROGRAMS = $(TESTS)

test_unparametric_field_SOURCES =		\
	test-unparametric-field.cpp		\
	test-generic.h			\
	test-common.C	test-common.h
	
test_ntl_RR_SOURCES =		\
	test-ntl-RR.cpp		\
	test-generic.h			\
	test-common.C	test-common.h
	
test_ntl_ZZ_p_SOURCES =		\
	test-ntl-ZZ_p.cpp		\
	test-generic.h			\
	test-common.C	test-common.h
	
test_ntl_zz_p_SOURCES =		\
	test-ntl-zz_p.cpp		\
	test-generic.h			\
	test-common.C	test-common.h
	
test_large_modular_SOURCES =		\
	test-large-modular.C		\
	test-generic.h			\
	test-common.C	test-common.h

test_vector_domain_SOURCES = 		\
	test-vector-domain.C		\
	test-common.C	test-common.h

test_dense_matrix_SOURCES = 		\
	test-dense-matrix.C		\
	test-common.C	test-common.h

test_sparse0_SOURCES = 		\
	test-sparse0.C		\
	test-common.C	test-common.h

test_diagonal_SOURCES = 		\
	test-diagonal.C			\
	test-generic.h			\
	test-common.C	test-common.h

test_minpoly_SOURCES = 			\
	test-minpoly.C			\
	test-common.C	test-common.h

test_det_SOURCES = 			\
	test-det.C			\
	test-common.C	test-common.h

test_rank_SOURCES =			\
	test-rank.C			\
	test-common.C	test-common.h

test_inverse_SOURCES =			\
	test-inverse.C			\
	test-common.C	test-common.h

test_unparametric_field_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
test_ntl_RR_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
test_ntl_ZZ_p_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
test_ntl_zz_p_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
test_large_modular_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
test_vector_domain_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
test_sparse0_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
test_dense_matrix_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
test_diagonal_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
test_minpoly_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
test_det_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
test_rank_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
test_inverse_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
#	test-rational-solver-adaptive.C \
#	test-common.h

#test_smith_form_adaptive_SOURCES=    \
#	test-common.h	\
#	test-smith-form-adaptive.C
