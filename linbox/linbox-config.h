/* Copyright (C) 2013 the members of the LinBox group
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *
 * This file is part of the LinBox library.
 *
 * ========LICENCE========
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * LinBox is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 */

#ifndef __LINBOX_linbox_config_H
#define __LINBOX_linbox_config_H

/** @file linbox-config.h
 * @brief linbox base configuration file
 */
#include <fflas-ffpack/fflas-ffpack-config.h>
#include "linbox/config.h"
// #include "linbox-configuration.h"

#ifndef INT32_MAX
#define INT32_MAX (2147483647L)
#endif

//#include <cstdint>
// #include <givaro/givconfig.h>
// #include <gmp++/gmp++.h>
#include <cfloat> // BB : needed on some rare platforms...
#include <cmath>

#include <cctype>
#include <iostream>

using std::ptrdiff_t;


#ifdef __LINBOX_HAVE_STDINT_H
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
// else ??
#endif
#include <stdint.h>
#ifndef INT32_MAX
#error "INT32_MAX is not defined. It should at least be defined in Givaro..."
#endif
#endif

/* Define if sse instructions are supported */
#ifdef __SSE__
#define __LINBOX_HAVE_SSE_INSTRUCTIONS  1
#endif

/* Define if sse2 instructions are supported */
#ifdef __SSE2__
#define __LINBOX_HAVE_SSE2_INSTRUCTIONS  1
#endif

/* Define if sse3 instructions are supported */
#ifdef __SSE3__
#define __LINBOX_HAVE_SSE3_INSTRUCTIONS  1
#endif

/* Define if sse4.1 instructions are supported */
#ifdef __SSE4_1__
#define __LINBOX_HAVE_SSE4_1_INSTRUCTIONS  1
#endif

/* Define if sse4.2 instructions are supported */
#ifdef __SSE4_2__
#define __LINBOX_HAVE_SSE4_2_INSTRUCTIONS  1
#endif


/* 256 SIMD registers are not supported by gcc on CYGWIN
 * See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=54412
 */
#if not defined(__CYGWIN__) or not defined(__GNUC__)

/* Define if avx instructions are supported */
#ifdef __AVX__
#define __LINBOX_HAVE_AVX_INSTRUCTIONS  1
#endif

/* Define if avx2 instructions are supported */
#ifdef __AVX2__
#define __LINBOX_HAVE_AVX2_INSTRUCTIONS  1
#endif

/* Define if avx512f instructions are supported */
#ifdef __AVX512F__
#define __LINBOX_HAVE_AVX512F_INSTRUCTIONS  1
#endif

/* Define if fma instructions are supported */
#ifdef __FMA__
#define __LINBOX_HAVE_FMA_INSTRUCTIONS  1
#endif

#endif // CYGWIN and GCC

namespace LinBox {

	typedef ptrdiff_t index_t;

	const int BlasBound = 1 << 26;

	//! used to separate BLAS2 and BLAS3 operations
	struct ContainerCategories {
		struct Any {} ;
		struct Vector : public Any {};
		struct Matrix : public Any {};
		struct Other  : public Any {};
	} ;

	//! Trait for the Category
	template<class Container>
	struct ContainerTraits {
		typedef ContainerCategories::Any  ContainerCategory ;
	} ;
}


#endif // __LINBOX_linbox_config_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
