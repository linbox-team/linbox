/* Copyright (C) LinBox
 * Written by bds and zw
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

#ifndef __LINBOX_smith_form_adaptive_H
#define __LINBOX_smith_form_adaptive_H

/*! @file algorithms/smith-form-adaptive.h
 * @ingroup algorithms
 * Implement the adaptive algorithm for Smith form computation
 */

#include <vector>
#include "linbox/integer.h"
#include "linbox/matrix/dense-matrix.h"

namespace LinBox
{

	class SmithFormAdaptive {
	public:

		static const long prime[];// = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};

		static const int NPrime;// = 25;

		/* Compute the local smith form at prime p, when modular (p^e) fits in long
		 * Should work with SparseMatrix and BlasMatrix
		 */
		template <class Matrix>
		static void compute_local_long (BlasVector<PID_integer>& s, const Matrix& A, long p, long e);

		/* Compute the local smith form at prime p, when modular (p^e) doesnot fit in long
		 * Should work with SparseMatrix and BlasMatrix
		 */
		template <class Matrix>
		static void compute_local_big (BlasVector<PID_integer>& s, const Matrix& A, long p, long e);

		/* Compute the local smith form at prime p
		*/
		template <class Matrix>
		static void compute_local (BlasVector<PID_integer>& s, const Matrix& A, long p, long e);

		/* Compute the k-smooth part of the invariant factor, where k = 100.
		 * @param sev is the exponent part ...
		 * By local smith form and rank computation
		 * Should work with SparseMatrix and BlasMatrix
		 */
		template <class Matrix>
		static void smithFormSmooth (BlasVector<PID_integer>& s, const Matrix& A, long r, const std::vector<long>& sev);

		/* Compute the k-rough part of the invariant factor, where k = 100.
		 * By EGV+ algorithm or Iliopoulos' algorithm for Smith form.
		 * Should work with BlasMatrix
		 */
		template <class Matrix>
		static void smithFormRough  (BlasVector<PID_integer>& s, const Matrix& A, integer m );

		/* Compute the Smith form via valence algorithms
		 * Compute the local Smith form at each possible prime
		 * r >= 2;
		 * Should work with SparseMatrix and BlasMatrix
		 */
		template <class Matrix>
		static void smithFormVal (BlasVector<PID_integer>&s, const Matrix& A, long r, const std::vector<long>& sev);

		/** \brief Smith form of a dense matrix by adaptive algorithm.
		 *
		 * Compute the largest invariant factor, then, based on that,
		 * compute the rough and smooth part, separately.
		 * Should work with SparseMatrix and BlasMatrix
		 */
		template <class Matrix>
		static void smithForm (BlasVector<PID_integer>& s, const Matrix& A);
		/** Specialization for dense case*/
		// template <class IRing>
		// static void smithForm (BlasVector<PID_integer>& s, const BlasMatrix<IRing>& A);
		template <class IRing, class _Rep>
		static void smithForm (BlasVector<PID_integer>& s, const BlasMatrix<IRing, _Rep>& A);

	};
	const long SmithFormAdaptive::prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
	const int SmithFormAdaptive::NPrime = 25;
}

#include "linbox/algorithms/smith-form-adaptive.inl"
#endif //__LINBOX_smith_form_adaptive_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
