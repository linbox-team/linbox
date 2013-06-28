/* linbox/algorithms/gauss-solve.inl
 * Copyright (C) LinBox 2008
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 * Time-stamp: <23 Mar 12 17:33:46 Jean-Guillaume.Dumas@imag.fr>
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

#ifndef __LINBOX_gauss_solve_INL
#define __LINBOX_gauss_solve_INL

// #include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/triangular-solve.h"
#include "linbox/blackbox/permutation.h"

namespace LinBox
{


	template <class _Field>
	template <class Matrix, class Perm, class Vector1, class Vector2> inline Vector1&
	GaussDomain<_Field>::solve(Vector1& x, Vector1& w, unsigned long Rank, const Perm& Q, const Matrix& L, const Matrix& U, const Perm& P, const Vector2& b)  const
	{

		Vector2 y(U.field(),U.rowdim()), v(U.field(),U.rowdim());

		Q.applyTranspose(y, b);

		lowerTriangularUnitarySolve(v, L, y);

		upperTriangularSolve(w, U, v);

		return P.applyTranspose(x, w);
	}

	template <class _Field>
	template <class Matrix, class Vector1, class Vector2> inline Vector1&
	GaussDomain<_Field>::solvein(Vector1& x, Matrix& A, const Vector2& b)  const
	{

		typename Field::Element Det;
		unsigned long Rank;
		Matrix L(field(), A.rowdim(), A.rowdim());
		Permutation<Field> Q((int)A.rowdim(),field());
		Permutation<Field> P((int)A.coldim(),field());

		this->QLUPin(Rank, Det, Q, L, A, P, A.rowdim(), A.coldim() );

		// Sets solution values to 0 for coldim()-Rank columns
		// Therefore, prune unnecessary elements
		// in those last columns of U
		for(typename Matrix::RowIterator row=A.rowBegin();
		    row != A.rowEnd(); ++row) {
			if (row->size()) {
				size_t ns=0;
				for(typename Matrix::Row::iterator it = row->begin();
				    it != row->end(); ++it, ++ns) {
					if (it->first >= Rank) {
						row->resize(ns);
						break;
					}
				}
			}
		}

		Vector1 w(A.field(),A.coldim());

		for(typename Vector1::iterator it=w.begin()+(ptrdiff_t)Rank;it!=w.end();++it)
			field().init(*it,0);

		return this->solve(x, w, Rank, Q, L, A, P, b);
	}

	template <class _Field>
	template <class Matrix, class Vector1, class Vector2, class Random> inline Vector1&
	GaussDomain<_Field>::solvein(Vector1& x, Matrix& A, const Vector2& b, Random& generator)  const
	{
		typename Field::Element Det;
		unsigned long Rank;
		Matrix L(field(), A.rowdim(), A.rowdim());
		Permutation<Field> Q((int)A.rowdim(),field());
		Permutation<Field> P((int)A.coldim(),field());

		this->QLUPin(Rank, Det, Q, L, A, P, A.rowdim(), A.coldim() );

                Vector1 w(A.field(),A.coldim());
		for(typename Vector1::iterator it=w.begin()+(ptrdiff_t)Rank;it!=w.end();++it)
			generator.random( *it );

		return this->solve(x, w, Rank, Q, L, A, P, b);
	}


} // namespace LinBox

#endif // __LINBOX_gauss_solve_INL


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

