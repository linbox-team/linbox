includedir=$(prefix)/include/linbox/matrix

BASIC_HDRS =			\
	archetype.h		\
	sparse.h		\
	sparse.inl		\
	transpose.h		\
	dense.h  		\
	dense.inl		\
	dense-submatrix.h	\
	dense-submatrix.inl	\
	dense-rows-matrix.h

NTL_HDRS =
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
	archetype.dxx		\
	matrices.dxx

