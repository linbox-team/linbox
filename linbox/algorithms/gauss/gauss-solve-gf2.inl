/* linbox/algorithms/gauss-solve-gf2.inl
 * Copyright (C) LinBox 2009
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 * Time-stamp: <13 Nov 17 16:59:23 Jean-Guillaume.Dumas@imag.fr>
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
	Vector1& GaussDomain<GF2>::solve(Vector1& x, Vector1& w, size_t Rank,
					 const Perm& Q, const SparseSeqMatrix& L,
					 const SparseSeqMatrix& U, const Perm& P,
					 const Vector2& b) const
	{

                const GF2 F2;
		Vector2 y(F2, U.rowdim()), v(F2, U.rowdim());

		Q.applyTranspose(y, b);

		lowerTriangularUnitarySolveBinary(v, L, y);

		upperTriangularSolveBinary(w, U, v);

		return P.applyTranspose(x, w);
	}

	template <class SparseSeqMatrix, class Vector1, class Vector2>
	Vector1& GaussDomain<GF2>::solveInPlace(Vector1& x,
					   SparseSeqMatrix        &A,
					   const Vector2& b) const
	{

		typename GF2::Element Det;
		size_t Rank;
		const GF2 F2;
		SparseSeqMatrix L(F2, A.rowdim(), A.rowdim());
		Permutation<GF2> Q(F2,(int)A.rowdim());
		Permutation<GF2> P(F2,(int)A.coldim());

		this->QLUPin(Rank, Det, Q, L, A, P, A.rowdim(), A.coldim() );

		Vector1 w(F2, A.coldim());

        for(typename Vector1::iterator it=w.begin()+(ptrdiff_t)Rank;it!=w.end();++it)
				F2.assign(*it,F2.zero);

		return this->solve(x, w, Rank, Q, L, A, P, b);
	}

	template <class SparseSeqMatrix, class Vector1, class Vector2, class Random>
	Vector1& GaussDomain<GF2>::solveInPlace(Vector1& x,
					   SparseSeqMatrix        &A,
					   const Vector2& b,
                       Random& generator) const
	{

		typename GF2::Element Det;
		size_t Rank;
		const GF2 F2;
		SparseSeqMatrix L(F2, A.rowdim(), A.rowdim());
		Permutation<GF2> Q((int)A.rowdim(),F2);
		Permutation<GF2> P((int)A.coldim(),F2);

		this->QLUPin(Rank, Det, Q, L, A, P, A.rowdim(), A.coldim() );

		Vector1 w(F2, A.coldim());

        for(typename Vector1::iterator it=w.begin()+Rank;it!=w.end();++it)
            generator.random( *it );

		return this->solve(x, w, Rank, Q, L, A, P, b);
	}

} // namespace LinBox

#endif // __LINBOX_gauss_solve_gf2_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
