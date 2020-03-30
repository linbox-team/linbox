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

AM_CPPFLAGS=-DLinBoxSrcOnly
AM_CXXFLAGS=$(LINBOX_CXXFLAGS) $(FFLAS_FFPACK_CFLAGS) $(LINBOX_DEPS_CFLAGS)
LDADD=$(LINBOX_DEPS_LIBS) $(FFLAS_FFPACK_LIBS) $(LINBOX_LDFLAGS)

pkgincludesubdir=$(pkgincludedir)/solutions

SUBDIRS=solve echelon

pkgincludesub_HEADERS=          \
    charpoly.h                  \
    det.h                       \
    echelon.h                   \
    getentry.h                  \
    getentry.inl                \
    is-positive-definite.h      \
    is-positive-semidefinite.h  \
    methods.h                   \
    minpoly.h                   \
    nullspace.h                 \
    rank.h                      \
    rank.inl                    \
    smith-form.h                \
    solve.h                     \
    trace.h                     \
    trace.inl                   \
    valence.h			\
    hadamard-bound.h            \
    constants.h	                \
    solution-tags.h

#    rankInPlace.h
#    wiedemann-rank.h

EXTRA_DIST =                \
    charpoly.doxy           \
    determinant.doxy        \
    minpoly.doxy            \
    nullspace.doxy          \
    rank.doxy               \
    reducedforms.doxy       \
    solutions.doxy          \
    systemsolving.doxy
