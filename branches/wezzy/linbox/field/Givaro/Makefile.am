# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
pkgincludesubdir=$(pkgincludedir)/field/Givaro

GIVARO_HDRS =           \
    givaro-zpz.h        \
    givaro-montg.h      \
    givaro-zpz.inl      \
    givaro-gfq.h        \
    givaro-extension.h  \
    givaro-field.h      \
    givaro-field-gf2.h  \
    givaro-rational.h

if LINBOX_HAVE_GIVARO
USE_GIVARO_HDRS = $(GIVARO_HDRS)
endif

pkgincludesub_HEADERS =     \
    $(USE_GIVARO_HDRS)

#  EXTRA_DIST =            \
    #  givaro.doxy
