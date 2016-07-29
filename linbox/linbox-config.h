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
#include "config.h"
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

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
#define __LINBOX_HAVE_SSE4_1_INSTRUCTIONS
#else
#define __LINBOX_NO_SIMD
#endif

#ifdef __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS
#define __LINBOX_HAVE_AVX_INSTRUCTIONS
#endif

#ifdef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
#define __LINBOX_HAVE_AVX2_INSTRUCTIONS
#endif

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
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

