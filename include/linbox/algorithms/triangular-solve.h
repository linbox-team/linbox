/* ===================================================================
 * Copyright(C) 2008 LinBox
 * Triangular Solve
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
 * Time-stamp: <16 Jun 10 14:21:18 Jean-Guillaume.Dumas@imag.fr>
 * ===================================================================
 */
#ifndef __LINBOX_triangular_solve_H
#define __LINBOX_triangular_solve_H

#include "linbox/vector/vector-domain.h"

namespace LinBox
{
	template <class _Matrix, class Vector1, class Vector2> Vector1&
	upperTriangularSolve (Vector1& x,
			      const _Matrix  &U,
			      const Vector2& b)
	{
		linbox_check( x.size() == U.coldim() );
		linbox_check( b.size() == U.rowdim() );
		typedef _Matrix Matrix;
		typedef typename Matrix::Field Field;
		const Field& F = U.field();

		commentator().start ("Sparse Elimination Upper Triangular Solve", "utrsm");

		typename Vector2::const_iterator vec=b.begin();
		typename Vector1::iterator res=x.begin();
		typename Matrix::ConstRowIterator row=U.rowBegin();

		// Assume U has form (U1, X | 0, 0), where U1 is invertible.
		// Discover the rank of U so as to use the bottom of x to seed the solution.
		for(row = U.rowEnd()-1; row >= U.rowBegin() && row->size() == 0; --row);
		++row; // now points to first zero row of U.
		res += row - U.rowBegin();
		vec += row - U.rowBegin();

		bool consistant = true;
		for(typename Vector2::const_iterator bcheck=vec; bcheck != b.end(); ++bcheck) {
			if( ! F.isZero(*bcheck) ) {
				consistant = false;
				break;
			}
		}
		if (consistant) {
			--vec; --res; --row;

			VectorDomain<Field> VD(F);
			for( ; row != U.rowBegin(); --row, --vec, --res) {
				F.assign(*res, F.zero);
				if (row->size()) {
					typename Field::Element tmp;
					VD.dot(tmp, *row, x);
					F.negin(tmp);
					F.addin(tmp,*vec);
					F.divin(tmp,row->front().second);
					F.assign(*res,tmp);
				}
				else {
					// Consistency check
					if( ! F.isZero(*vec) ) {
						consistant = false;
						break;
					}
				}
			}

			F.assign(*res, F.zero);
			if (row->size()) {
				typename Field::Element tmp;
				VD.dot(tmp, *row, x);
				F.negin(tmp);
				F.addin(tmp,*vec);
				F.divin(tmp,row->front().second);
				F.assign(*res,tmp);
			}
			else {
				// Consistency check
				if( ! F.isZero(*vec) ) consistant = false;
			}
		}
		if (! consistant) throw LinboxError ("upperTriangularSolve returned INCONSISTENT");

		commentator().stop ("done", NULL, "utrsm");
		return x;
	}

	// Suppose x and b are vectors of pairs <index,value>
	// first rank rows of U are upper triangular and full rank
	template <class _Matrix, class Vector1, class Vector2> Vector1&
	upperTriangularSparseSolve (Vector1& x,
				    unsigned long rank,
				    const _Matrix  &U,
				    const Vector2& b)
	{
		commentator().start ("SparseElim UpperTriang Sparse Solve", "uSPt");

		x.resize(0);

		if (b.size() != 0) {

			typedef _Matrix Matrix;
			typedef typename Matrix::Field Field;
			const Field& F = U.field();



			typename Vector2::const_iterator vec=b.begin();
			vec += b.size(); --vec;

			typename Matrix::ConstRowIterator row=U.rowBegin();
			row += (ptrdiff_t)rank;
		       	--row;

			VectorDomain<Field> VD(F);


			long i=(long)rank;
			for(--i; (vec >= b.begin()) && (i>=0); --i,--row) {
				if (row->size()) {
					typename Field::Element tmp;
					VD.dot(tmp, *row, x); // x is sparse also
					F.negin(tmp);
					if (static_cast<long>(vec->first) == i) {
						F.addin(tmp,vec->second);
						--vec;
					}
					if (! F.isZero(tmp)) {
						F.divin(tmp,row->front().second);
						x.insert(x.begin(), typename Vector1::value_type((unsigned)i, tmp));
					}
				}
			}
			for(; i>=0; --i,--row) {
				if (row->size()) {
					typename Field::Element tmp;
					VD.dot(tmp, *row, x);
					if (! F.isZero(tmp)) {
						F.negin(tmp);
						F.divin(tmp,row->front().second);
						x.insert(x.begin(), typename Vector1::value_type((unsigned)i, tmp));
					}
				}
			}
		}
		//         if (! consistant) throw LinboxError ("upperTriangularSparseSolve returned INCONSISTENT");

		commentator().stop ("done", NULL, "uSPt");

		return x;
	}


	template <class _Matrix, class Vector1, class Vector2> Vector1&
	lowerTriangularUnitarySolve (Vector1& x,
				     const _Matrix  &L,
				     const Vector2& b)
	{
		linbox_check( b.size() == L.coldim() );
		typedef _Matrix Matrix;
		typedef typename Matrix::Field Field;
		const Field& F = L.field();

		commentator().start ("Sparse Elimination Lower Triangular Unitary Solve", "ltrsm");

		typename Vector2::const_iterator vec=b.begin();
		typename Vector1::iterator res=x.begin();
		typename Matrix::ConstRowIterator row=L.rowBegin();

		VectorDomain<Field> VD(F);
		for( ; row != L.rowEnd(); ++row, ++vec, ++res) {
			F.assign(*res, F.zero);
			typename Field::Element tmp;
			VD.dot(tmp, *row, x);
			F.negin(tmp);
			F.addin(tmp,*vec);
			F.assign(*res,tmp);
		}

		commentator().stop ("done", NULL, "ltrsm");
		return x;
	}
}
#endif //__LINBOX_triangular_solve_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

