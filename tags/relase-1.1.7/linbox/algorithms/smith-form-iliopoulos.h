/* Copyright (C) 2010 LinBox
 *
 *  Author: Zhendong Wan
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

#ifndef __LINBOX_smith_form_iliopoulos_H
#define __LINBOX_smith_form_iliopoulos_H

#include <linbox/util/debug.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/submatrix-traits.h>

namespace LinBox 
{
	
/** \brief This is Iliopoulos' algorithm do diagonalize.

 * Compute Smith Form by elimination modulo m, for some modulus m such 
 * as S(n), the last invariant factor.  
 * The elimination method is originally described in 
 * "Worst Case Complexity Bounds on Algorithms for computing the Canonical 
 *  Structure of Finite Abelian Groups and the Hermite and Smith Normal 
 * Forms of an Integer Matrix", by Costas Iliopoulos.
*/

	class SmithFormIliopoulos{
		

	protected:
		/** \brief eliminationRow will make the first row (*, 0, ..., 0)
		 *  by col operations.
		 *  It is the implementation of Iliopoulos algorithm
		 */
		template<class Matrix, class Ring>
		static Matrix& eliminationRow (Matrix& A, const Ring& r) {
						      
			
			if (A. coldim() <= 1 || A. coldim()  == 0) return A;
			

			//typedef typename Matrix::Field Field;
			
			typedef typename Ring::Element Element;

		
			VectorDomain<Ring> vd (r);

			// some tempory variable
			typename Matrix::RowIterator cur_r, tmp_r;
			
			typename Matrix::ColIterator cur_c, tmp_c;
			typename Matrix::Row::iterator row_p1, row_p2;
			//typename Matrix::Col::iterator col_p1, col_p2;
			
			cur_c = A. colBegin();

			cur_r = A. rowBegin();
			
			row_p1 = cur_r -> begin();
			
			// if A[0][0] is coprime to d
			if (r. isUnit( *row_p1)) {

				if (! r. isOne (* row_p1)) {
					
					Element s;

					r. inv (s, *row_p1);

					vd. mulin(*cur_c, s);

				}
				
			}
			
			// A[0][0] is not a unit
			else {
					
				// make A[0][0] = 0
				if ( !r.isZero(*row_p1)) {
					
					row_p2 = row_p1 + 1;
					
					Element y1, y2;

					Element g, s, t;
					
					r. dxgcd (g, s, t, y2, y1, *row_p1, *row_p2);

					r. negin (y1);
					
					tmp_c = cur_c + 1;
					
					std::vector<Element> tmp1 (A.rowdim()), tmp2 (A.rowdim());

					vd. mul (tmp1, *cur_c, y1);

					vd. axpyin (tmp1, y2, *tmp_c);

					vd. mul (tmp2, *cur_c, s);

					vd. axpyin (tmp2, t, *tmp_c);

					vd. copy (*cur_c,  tmp1);

					vd. copy (*tmp_c, tmp2);

					if (!r.isZero (*(cur_c -> begin()))) {

						Element q;

						r. div (q, *(cur_c -> begin()), g);

						r. negin (q);
							  
						vd. axpyin (*cur_c, q, *tmp_c);
					}
				}
				
				// matrix index is 0-based
				std::vector<Element> tmp_v(A.coldim());
				
				typename std::vector<Element>::iterator p1, p2;

				r. init(tmp_v[0], 1);

				p1 = tmp_v.begin() + 1;
				p2 = tmp_v.begin() + 1;
				
				row_p2 = row_p1 + 1; 
				
				Element g, s;
				
				r.assign(g, *row_p2); ++ row_p2;
				
				r.init(*p1, 1); ++ p1;
				
				for (; row_p2 != cur_r -> end(); ++ row_p2, ++ p1) {
					
					r.xgcd(g, s, *p1, g, *row_p2);
					
					if (!r.isOne(s))
						
						for (p2 = tmp_v.begin() + 1; p2 != p1; ++ p2) 
							
							r. mulin (*p2, s);
					
						
					
				}
				
				// no pivot found
				if (r.isZero(g)) return A;

				
				for (tmp_r = cur_r; tmp_r != A.rowEnd(); ++ tmp_r) 
					
					vd. dot (*(tmp_r -> begin()), *tmp_r, tmp_v);
					
				
			}
				
			
			// after finding the pivot
			// column operation to make A[p][j] = 0, where k < j


			Element g, tmp;
			
			r. assign (g, *(cur_c -> begin()));

			for (tmp_c = cur_c + 1; tmp_c != A.colEnd(); ++ tmp_c) {
							       
				// test if needing to update 
				if (!r. isZero (*(tmp_c -> begin()))) {
					
					r.div (tmp, *(tmp_c -> begin()), g);
					
					r.negin(tmp);
					
					vd. axpyin (*tmp_c, tmp, *cur_c);

				}
			}


			return A;
		}
		

		
		/** \brief eliminationCol will make the first col (*, 0, ..., 0)
		 *  by elementary row operation.
		 *  It is the implementation of Iliopoulos algorithm
		 */
		template<class Matrix, class Ring>
		static Matrix& eliminationCol (Matrix& A, const Ring& r) {

			if((A.rowdim() <= 1) || (A.rowdim() == 0)) return A;
			
			//typedef typename Matrix::Field Field;
			typedef typename Ring::Element Element;
			
			//Field r (A.field());

			VectorDomain<Ring> vd(r);
			
			typename Matrix::ColIterator cur_c, tmp_c;
			typename Matrix::RowIterator cur_r, tmp_r;
			
			//typename Matrix::Row::iterator row_p1, row_p2;
			typename Matrix::Col::iterator col_p1, col_p2;
						
			
			cur_c = A.colBegin();
			cur_r = A.rowBegin();

			col_p1 = cur_c -> begin();


			// If A[0][0] is a unit
			if (r.isUnit (*col_p1) ) {

				if (! r. isOne ( *col_p1)) {
				
					
					Element s;

					r. inv (s, *col_p1);
					
					vd. mulin (*cur_r, s);
					
				}
				
			}
			else {
				// Make A[0][0] = 0;
				if (!r.isZero(*col_p1)) {
					
					Element g, s, t, y1, y2;
					
					std::vector<Element> tmp1(A.coldim()), tmp2(A.coldim());
					
					col_p2 = col_p1 + 1;
					
					r.dxgcd(g, s, t, y2, y1, *col_p1, *col_p2);

					r. negin (y1);
					

					tmp_r = cur_r + 1;
					
					vd. mul (tmp1, *cur_r, y1);
					
					vd. axpyin (tmp1, y2, *tmp_r);

					vd. mul (tmp2, *cur_r, s);

					vd. axpyin (tmp2, t, *tmp_r);

					vd. copy (*cur_r, tmp1);

					vd. copy (*tmp_r, tmp2);

					if ( !r. isZero (* (cur_r ->begin() ) ) ) {

						Element q;

						r. div (q, *(cur_r -> begin() ), g);

						r. negin (q);

						vd. axpyin (*cur_r, q, *tmp_r);
					}

				}


				// matrix index is 0-based
				std::vector<Element> tmp_v (A.rowdim());
					
				typename std::vector<Element>::iterator p1, p2;
				
				Element g, s;
				
				col_p2 = col_p1 + 1;
				
				r.assign (g, *col_p2); ++ col_p2;
				
				r. init (tmp_v[0], 1);
				
				p1 = tmp_v.begin() + 1;
				
				r.init(*p1,1); ++ p1;

				for(; col_p2 != cur_c -> end(); ++ col_p2, ++ p1) {
					
					r.xgcd (g, s, *p1, g, *col_p2);
					
					if (! r.isOne(s)) 
						for (p2 = tmp_v.begin() + 1; p2 != p1; ++ p2) 
							
							r. mulin (*p2, s);

					
				}
				
				if (r.isZero(g))  return A;

				// no pivot found
				for (tmp_c = cur_c; tmp_c != A.colEnd(); ++ tmp_c) 
					
					vd. dot ( *(tmp_c -> begin()), *tmp_c, tmp_v);
				
			}			


			// A pivot is found

			Element g, tmp;

			r. assign (g, *( cur_r -> begin()));
			
			for (tmp_r = cur_r + 1; tmp_r != A.rowEnd(); ++ tmp_r) {
					
				if (! r.isZero(*(tmp_r -> begin() ) ) ) {
					
					r.div (tmp, *(tmp_r -> begin()), g);
					
					r.negin (tmp);
					
					vd. axpyin (*tmp_r, tmp, *cur_r);

				}
			}


			return A;

		}

		template<class Matrix, class Ring>
		static bool check(const Matrix& A, const Ring& r) {
				
			
			//typedef typename Matrix::Ring Field;
			typedef typename Ring::Element Element;
			
			typename Matrix::ConstRowIterator cur_r;
			typename Matrix::ConstRow::const_iterator row_p;
			
			Element tmp;
			
			cur_r = A.rowBegin();
			row_p = cur_r -> begin();
				
			tmp = * (A.rowBegin() -> begin());
			
			if (r.isZero(tmp)) return true;

			for (++ row_p; row_p != cur_r -> end(); ++ row_p ) {
				
				if (!r. isDivisor (tmp, *row_p))

					return false;
			}
			
			return true;
		}
				
		/** \brief Diagonalize the matrix A.
		 */
		template<class Matrix, class Ring>
		static Matrix& diagonalizationIn(Matrix& A, const Ring& r) {
			
			if (A.rowdim() == 0 || A.coldim() == 0) return A;

			
			do {
				
				eliminationRow (A, r);

				eliminationCol (A, r);
			}

			while (!check(A, r));

			typename SubMatrixTraits<Matrix>::value_type 
				sub(A, (unsigned int)1, (unsigned int)1, 
				    A.rowdim() - 1, A.coldim() - 1);

			diagonalizationIn(sub, r);

			return A;
		}


	public:

		template<class Matrix>
		static  Matrix& smithFormIn(Matrix& A) {
			
			typedef typename Matrix::Field Ring;
			typedef typename Ring::Element Element;
			
			Ring r (A.field());

			typename Matrix::RowIterator row_p;
			
			Element zero, one;

			r. init (zero, 0);

			r. init (one, 1);
			
			diagonalizationIn(A, r);

			int min = A.rowdim() <= A.coldim() ? A.rowdim() : A.coldim();			

			int i, j;
			
			Element g;
			
			for (i = 0; i < min; ++ i) {
				
				for ( j = i + 1; j < min; ++ j) {
					
					if (r. isUnit(A[i][i]))  break;
						
					else if (r. isZero (A[j][j])) continue;
					
					else if (r. isZero (A[i][i])) {
						std::swap (A[i][i], A[j][j]);
					}
					
					else {
						r. gcd (g, A[j][j], A[i][i]);
						
						r. divin (A[j][j], g);
						
						r. mulin (A[j][j], A[i][i]);

						r. assign (A[i][i], g);
					}
				}
				r. normalIn (A[i][i]);
			}
			
			return A;
			
		}
		
	};
	
	
}
				
#endif //__LINBOX_smith_form_iliopoulos_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
