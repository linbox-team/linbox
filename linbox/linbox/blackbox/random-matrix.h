/* Copyright (C) 2010 LinBox
 *  Author: Zhendong Wan
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */



#ifndef __LINBOX_random_matrix_H
#define __LINBOX_random_matrix_H

#include <linbox/blackbox/blackbox-interface.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/sparse.h>

namespace LinBox 
{

	class RandomMatrix : public  BlackboxInterface 
	{

	public:

		/** generates random matrices used in EGV and EGV+ algorithm

		 * [I, R] or [I, R]^t, where R is a random matrix.
		 General case.
		 */
		template <class Blackbox, class Field>
		static Blackbox*& randomMatrix (Blackbox* &, const Field& f, 
						int rowdim, int coldim);


		// constructor a random dense matrix, whose entries are random
		template<class Field>
		static DenseMatrix<Field>*& randomMatrix( DenseMatrix<Field>*& Ap,
							  const Field& f,
							  int rowdim, int coldim ) {

			Ap = new DenseMatrix<Field>(f, rowdim, coldim);
			typename DenseMatrix<Field>::RawIterator Ap_p;
			typename Field::Element zero, one, elt;
			f. init (one, 1); f. init (zero, 0);

			for (Ap_p = Ap -> rawBegin(); Ap_p != Ap -> rawEnd(); ++ Ap_p)
				f. assign (*Ap_p, zero);

			if (rowdim < coldim) 
				for (int i = 0; i < rowdim; ++ i) {
					Ap -> setEntry (i, i, one);
					for (int j = rowdim; j < coldim; ++ j){
						f. init (elt, rand()%10);
						Ap -> setEntry (i, j, elt);
					}
				}
			else 
				for (int i = 0; i < coldim; ++ i) {
					Ap -> setEntry (i, i, one);
					for (int j = coldim; j < rowdim; ++ j) {
						f. init (elt, rand()%10);
						Ap -> setEntry (j, i, elt);
					}
				}


			return Ap;
		}

		// constructor a very special random sparse matrix
		// [I, R] or [I, R}^t, where R is a sparse random matrix.
		template<class Field>
		static SparseMatrix<Field>*& randomMatrix( SparseMatrix<Field>*& Ap, 
							   const Field& f, 
							   int rowdim, int coldim) {

			Ap = new SparseMatrix<Field>(f, rowdim, coldim);

			const int m = rowdim < coldim ? rowdim : coldim;

			int i, j, k;

			typename Field::Element elt;

			f. init (elt, 1);

			for ( i = 0; i < m; ++ i) 

				Ap -> setEntry (i, i, elt);


			if ( m < rowdim ) {

				const int repeat = (rowdim - m) < 10 ? rowdim - m : 10;

				for ( i = m; i < rowdim; ++ i) {

					for ( k = 0; k < repeat; ++ k) {

						j = rand() % coldim;

						f.init(elt, rand() % 10 + 1);

						Ap -> setEntry (i, j, elt);
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

						Ap -> setEntry (i, j, elt);
					}
				}
			}

			else {}

			return Ap;

		}

		template<typename _Tp1> 
		struct rebind 
		{ typedef RandomMatrix other; };

	};
}


#endif //__LINBOX_random_matrix_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
