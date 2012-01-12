/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/gauss-solve.inl
 * Copyright (C) LinBox 2008
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 * Time-stamp: <21 Jan 10 17:01:22 Jean-Guillaume.Dumas@imag.fr>
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

#include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/triangular-solve.h"
#include "linbox/blackbox/permutation.h"

namespace LinBox
{


	template <class _Field>
	template <class Matrix, class Perm, class Vector1, class Vector2> inline Vector1&
	GaussDomain<_Field>::solve(Vector1& x, Vector1& w, unsigned long Rank, const Perm& Q, const Matrix& L, const Matrix& U, const Perm& P, const Vector2& b)  const
	{

		Vector2 y(U.rowdim()), v(U.rowdim());

		Q.applyTranspose(y, b);

		lowerTriangularUnitarySolve(v, L, y);

		upperTriangularSolve(w, U, v);

		return P.applyTranspose(x, w);
	}

	template <class _Field>
	template <class Matrix, class Perm, class Vector1, class Vector2> inline Vector1&
	GaussDomain<_Field>::solve(Vector1& x, unsigned long Rank, const Perm& Q, const Matrix& L, const Matrix& U, const Perm& P, const Vector2& b, bool randomsol)  const
	{

		Vector1 w(U.coldim());
		if (randomsol) {
			// Random solution is in output
			typename _Field::RandIter generator(_field);
			for(typename Vector1::iterator it=w.begin()+Rank;it!=w.end();++it)
				generator.random( *it );
		}
		else {
			for(typename Vector1::iterator it=w.begin()+Rank;it!=w.end();++it)
				_field.init(*it,0);
		}
		return this->solve(x, w, Rank, Q, L, U, P, b);
	}

	template <class _Field>
	template <class Matrix, class Vector1, class Vector2> inline Vector1&
	GaussDomain<_Field>::solvein(Vector1& x, Matrix& A, const Vector2& b, bool randomsol)  const
	{

		typename Field::Element Det;
		unsigned long Rank;
		Matrix L(_field, A.rowdim(), A.rowdim());
		Permutation<Field> Q((int)A.rowdim(),_field);
		Permutation<Field> P((int)A.coldim(),_field);

		this->QLUPin(Rank, Det, Q, L, A, P, A.rowdim(), A.coldim() );

		if (! randomsol) {
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
		}

		return this->solve(x, Rank, Q, L, A, P, b, randomsol);
	}

} // namespace LinBox

#endif // __LINBOX_gauss_solve_INL

