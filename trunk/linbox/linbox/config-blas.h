/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* config-blas.h
 * Copyright (C) 2005  Pascal Giorgi
 *               2007  Clement Pernet
 * Written by Pascal Giorgi <pgiorgi@uwaterloo.ca>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


#ifndef __LINBOX_config_blas_H
#define __LINBOX_config_blas_H

#include <fflas-ffpack/config-blas.h>
#include <fflas-ffpack/fflas-ffpack-config.h>

#ifdef __FFLASFFPACK_HAVE_BLAS
#define __LINBOX_HAVE_BLAS 1
#endif

#ifdef __FFLASFFPACK_HAVE_CBLAS
#define __LINBOX_HAVE_CBLAS 1
#endif

#ifdef __FFLASFFPACK_HAVE_LAPACK
#define __LINBOX_HAVE_LAPACK 1
#endif

#ifdef __FFLASFFPACK_HAVE_CLAPACK
#define __LINBOX_HAVE_CLAPACK 1
#endif

#endif //__LINBOX_config_blas_H
