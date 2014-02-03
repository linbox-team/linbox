/* linbox/algorithms/gauss-solve-gf2.inl
 * Copyright (C) LinBox 2009
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 * Time-stamp: <23 Mar 12 17:32:19 Jean-Guillaume.Dumas@imag.fr>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_gauss_solve_gf2_INL
#define __LINBOX_gauss_solve_gf2_INL

#include "linbox/algorithms/gauss-gf2.h"
#include "linbox/algorithms/triangular-solve-gf2.h"
#include "linbox/blackbox/permutation.h"

namespace LinBox
{


	template <class SparseSeqMatrix, class Perm, class Vector1, class Vector2>
	Vector1& GaussDomain<GF2>::solve(Vector1& x, Vector1& w, unsigned long Rank,
					 const Perm& Q, const SparseSeqMatrix& L,
					 const SparseSeqMatrix& U, const Perm& P,
					 const Vector2& b) const
	{

		Vector2 y(U.rowdim()), v(U.rowdim());

		Q.applyTranspose(y, b);

		lowerTriangularUnitarySolveBinary(v, L, y);

		upperTriangularSolveBinary(w, U, v);

		return P.applyTranspose(x, w);
	}

	template <class SparseSeqMatrix, class Vector1, class Vector2>
	Vector1& GaussDomain<GF2>::solvein(Vector1& x,
					   SparseSeqMatrix        &A,
					   const Vector2& b) const
	{

		typename GF2::Element Det;
		unsigned long Rank;
		const GF2 F2;
		SparseSeqMatrix L(F2, A.rowdim(), A.rowdim());
		Permutation<GF2> Q((int)A.rowdim(),F2);
		Permutation<GF2> P((int)A.coldim(),F2);

		this->QLUPin(Rank, Det, Q, L, A, P, A.rowdim(), A.coldim() );

		Vector1 w(A.coldim());

        for(typename Vector1::iterator it=w.begin()+(ptrdiff_t)Rank;it!=w.end();++it)
				F2.assign(*it,F2.zero);

		return this->solve(x, w, Rank, Q, L, A, P, b);
	}

	template <class SparseSeqMatrix, class Vector1, class Vector2, class Random>
	Vector1& GaussDomain<GF2>::solvein(Vector1& x,
					   SparseSeqMatrix        &A,
					   const Vector2& b,
                       Random& generator) const
	{

		typename GF2::Element Det;
		unsigned long Rank;
		const GF2 F2;
		SparseSeqMatrix L(F2, A.rowdim(), A.rowdim());
		Permutation<GF2> Q((int)A.rowdim(),F2);
		Permutation<GF2> P((int)A.coldim(),F2);

		this->QLUPin(Rank, Det, Q, L, A, P, A.rowdim(), A.coldim() );

		Vector1 w(A.coldim());

        for(typename Vector1::iterator it=w.begin()+Rank;it!=w.end();++it)
            generator.random( *it );

		return this->solve(x, w, Rank, Q, L, A, P, b);
	}

} // namespace LinBox

#endif // __LINBOX_gauss_solve_gf2_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

