# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
pkgincludesubdir=$(pkgincludedir)/NTL


NTL_HDRS =              \
    ntl-RR.h            \
    ntl-ZZ_p.h          \
    ntl-ZZ_pE.h         \
    ntl-ZZ_pX.h         \
    ntl-lzz_p.h         \
    ntl-ZZ.h            \
    ntl-lzz_pE.h        \
    ntl-lzz_pX.h        \
    ntl-GF2E.h


if LINBOX_HAVE_NTL
USE_NTL_HDRS = $(NTL_HDRS)
endif


pkgincludesub_HEADERS =     \
    $(USE_NTL_HDRS)

#  EXTRA_DIST =            \
    #  ntl.doxy
