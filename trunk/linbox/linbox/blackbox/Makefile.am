includedir=$(prefix)/include/linbox/blackbox

BASIC_HDRS =			\
	archetype.h		\
	butterfly.h		\
	companion.h		\
	compose.h		\
	dense.h  		\
	dense.inl		\
	diagonal.h		\
	dif.h			\
	direct-sum.h		\
	factory.h		\
	frobenius.h		\
	hilbert.h		\
	inverse.h		\
	matrix-blackbox.h	\
	moore-penrose.h		\
	null-matrix.h		\
	permutation.h		\
	scalar-matrix.h		\
	sparse.h		\
	sparse.inl		\
	sparse1.h		\
	submatrix.h		\
	sum.h			\
	transpose.h		\
	triplesbb.h		\
	zero-one.h		\
	zero-one.inl		



NTL_HDRS =			\
	ntl-hankel.h            \
	ntl-hankel.inl		\
	ntl-sylvester.h         \
	ntl-sylvester.inl	\
	ntl-toeplitz.h          \
	ntl-toeplitz.inl

GIVARO_HDRS = 

LIDIA_HDRS = 

if HAVE_NTL
USE_NTL_HDRS = $(NTL_HDRS)
endif

if HAVE_GIVARO
USE_GIVARO_HDRS = $(GIVARO_HDRS)
endif

if HAVE_LIDIA
USE_LIDIA_HDRS = $(LIDIA_HDRS)
endif

include_HEADERS =		\
	$(BASIC_HDRS)		\
	$(USE_NTL_HDRS)		\
	$(USE_GIVARO_HDRS)	\
	$(USE_LIDIA_HDRS)

EXTRA_DIST = 
	$(NTL_HDRS)		\
	$(GIVARO_HDRS)		\
	$(LIDIA_HDRS)		\
	blackbox.dxx		

