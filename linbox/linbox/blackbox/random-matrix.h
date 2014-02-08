/* Copyright (C) 2010 LinBox
 *  Author: Zhendong Wan
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



#ifndef __LINBOX_blackbox_random_matrix_H
#define __LINBOX_blackbox_random_matrix_H

#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/matrix/sparse-matrix.h"

namespace LinBox
{

	class RandomMatrix : public  BlackboxInterface {

	public:

		/** Generates random matrices used in EGV and EGV+ algorithm.

		 * [I, R] or [I, R]^t, where R is a random matrix.
		 * General case.
		 */
		template <class Blackbox, class Field>
		static Blackbox*& randomMatrix (Blackbox* &, const Field& f,
						int rowdim, int coldim);


		// constructor a random dense matrix, whose entries are random
		template<class Field>
		static BlasMatrix<Field>*& randomMatrix( BlasMatrix<Field>*& Ap,
							  const Field& f,
							  int rowdim, int coldim )
		{

			Ap = new BlasMatrix<Field>(f, rowdim, coldim);
			typename BlasMatrix<Field>::Iterator Ap_p;
			typename Field::Element  elt;

			for (Ap_p = Ap -> Begin(); Ap_p != Ap -> End(); ++ Ap_p)
				f. assign (*Ap_p, f.zero);

			if (rowdim < coldim)
				for (int i = 0; i < rowdim; ++ i) {
					Ap -> setEntry ((size_t)i,(size_t) i, f.one);
					for (int j = rowdim; j < coldim; ++ j){
						f. init (elt, rand()%10);
						Ap -> setEntry ((size_t)i, (size_t)j, elt);
					}
				}
			else
				for (int i = 0; i < coldim; ++ i) {
					Ap -> setEntry ((size_t)i,(size_t) i, f.one);
					for (int j = coldim; j < rowdim; ++ j) {
						f. init (elt, rand()%10);
						Ap -> setEntry ((size_t)j,(size_t) i, elt);
					}
				}


			return Ap;
		}


		// constructor a very special random sparse matrix
		// [I, R] or [I, R]^t, where R is a sparse random matrix.
		template<class Field>
		static SparseMatrix<Field>*& randomMatrix( SparseMatrix<Field>*& Ap,
							   const Field& f,
							   int rowdim, int coldim)
		{

			Ap = new SparseMatrix<Field>(f, rowdim, coldim);

			const int m = rowdim < coldim ? rowdim : coldim;

			int i, j, k;

			typename Field::Element elt;

			f. assign(elt, f.one);

			for ( i = 0; i < m; ++ i)

				Ap -> setEntry ((size_t)i, (size_t)i, elt);


			if ( m < rowdim ) {

				const int repeat = (rowdim - m) < 10 ? rowdim - m : 10;

				for ( i = m; i < rowdim; ++ i) {

					for ( k = 0; k < repeat; ++ k) {

						j = rand() % coldim;

						f.init(elt, rand() % 10 + 1);

						Ap -> setEntry ((size_t)i, (size_t)j, elt);
					}
				}
			}

			else if ( m < coldim ) {

				int offset = coldim - m;

				const int repeat = offset < 10 ? offset : 10;

				for ( i = 0; i < rowdim; ++ i){

					for ( k = 0; k < repeat; ++ k) {

						j = rand() % offset + m;

						f.init(elt, rand() % 10 + 1);

						Ap -> setEntry ((size_t)i, (size_t)j, elt);
					}
				}
			}

			else {}

			return Ap;

		}

		template<typename _Tp1>
		struct rebind {
			typedef RandomMatrix other;
		};

	};
}


#endif //__LINBOX_blabbox_random_matrix_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

