includedir=$(prefix)/include/linbox/blackbox

BASIC_HDRS =			\
	archetype.h		\
	factory.h		\
	butterfly.h		\
	hilbert.h		\
	compose.h		\
	permutation.h		\
	sparse1.h		\
	sparse.h		\
	sparse.inl		\
	sparse-base.h		\
	sparse-base.inl		\
	scalar-matrix.h		\
	submatrix.h		\
	inverse.h		\
	transpose.h		\
	sum.h			\
	diagonal.h		\
	moore-penrose.h		\
	dense.h  		\
	dense.inl		\
        nag-sparse.h            

NTL_HDRS =			\
	ntl-toeplitz.h

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
	blackboxes.dxx		\
	archetype.dxx		\
	auxillary.dxx

