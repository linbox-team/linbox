/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* Copyright (C) 2011 LinBox
 * Written by BB <brice.boyer@imag.fr>
 *
 *
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

/*! @file util/iml_wrapper.h
 * @ingroup util
 * @brief Wraper for IML routines and puts them in a dedicated namespace.
 * This wrapper includes iml.h definitions in an IML namespace.  The \#defines
 * found elsewhere in the library (eg X*ALLOC) are rewritten.  The functions
 * that are provided in the .a but not in the header are provided <code>extern
 * "C"</code>.
 */

#ifndef __LINBOX_util_iml_wrapper_H
#define __LINBOX_util_iml_wrapper_H

#ifdef __LINBOX_HAVE_IML
/*! @brief Namespace for IML routines.
 * This namespace is foreign to LinBox'.
 */
namespace IML {
	extern "C" {
#include "iml.h"

		FiniteField **
		findRNS (const FiniteField RNS_bound, const mpz_t mp_maxInter, long *length) ;
		FiniteField *
		repBound (const long len, const FiniteField *basis, const FiniteField *cmbasis) ;
		void
		basisProd (const long len, const FiniteField *basis, mpz_t mp_prod) ;
		void
		ChineseRemainder (const long len, const mpz_t mp_prod, \
				  const FiniteField *basis, const FiniteField *cmbasis, \
				  const FiniteField *bdcoeff, Double *Ac, mpz_t mp_Ac);
		void
		ChineseRemainderPos (const long len, const FiniteField *basis, \
				     const FiniteField *cmbasis, Double *Ac, mpz_t mp_Ac);

	}
	// #include <stdlib.h>
#ifndef XMALLOC
#define IML_XMALLOC(type, num)                                  \
	((type *) malloc ((num) * sizeof(type)))
#else
#define IML_XMALLOC(type, num)                                  \
	XMALLOC(type, num)
#endif

#ifndef XFREE
#define IML_XFREE(stale)                            do {        \
	if (stale) { free (stale);  stale = 0; }                    \
        } while (0)
#else
#define IML_XFREE(stale)                                        \
	XFREE(stale)
#endif

	extern "C" {
		void *xcalloc  (size_t num, size_t size);
		void *xmalloc  (size_t num);
		void *xrealloc (void *p, size_t num);

	}

}
#else
#error "you are using IML wrapper without IML available in LinBox."
#endif

#endif // __LINBOX_util_iml_wrapper_H
