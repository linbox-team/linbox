/* Copyright (C) 2003 LinBox
 * Written by Zhendong Wan
 *
 * modified by Pascal Giorgi 1/07/04
 * put the Field as template parameter 
 * and add Field F as a parameter
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


#ifndef __LINBOX_matrix_inverse_H
#define __LINBOX_matrix_inverse_H

#include <linbox/util/debug.h>
#include <linbox/util/error.h>
#include <vector>
#include <algorithm>

namespace LinBox 
{
	class MatrixInverse {
	
	public:

		/** \brief compute the inverse of a dense matrix, by Gaussian elimination.
		 *  The matrix should support ColIterator and RowIterator.
		 *  It returns 0, if an inverse is found, and  
		 *  returns 1, otherwise.
		 */
		template<class Field, class DenseMatrix>
		static  long matrixInverseIn(const Field& F, DenseMatrix& A) {
		    
			// check if A is a square matrix
			linbox_check(A.rowdim() == A. coldim());
	
			// PG 1/07/04 
			//typedef typename DenseMatrix::Field Field;

			// step1 PLU Inplcae, actually, LPA = U.
			std::vector<std::pair<int,int> > P;
			P.reserve (A.rowdim());

			typename DenseMatrix::RowIterator cur_r, tmp_r;
			typename DenseMatrix::ColIterator cur_c, tmp_c;
			typename DenseMatrix::Row::iterator cur_rp, tmp_rp;
			typename DenseMatrix::Col::iterator cur_cp, tmp_cp;
	
			std::vector<typename Field::Element> tmp_v (A.rowdim());

			typename Field::Element tmp_e;
			
			// PG 1/07/04
			//const Field F = A. field();
       
			int offset = 0;

			cur_r = A. rowBegin(); 
			cur_c = A. colBegin();
			for( ; cur_r != A. rowEnd(); ++ cur_r, ++ cur_c, ++ offset) { 
			//for(cur_r = A. rowBegin(), cur_c = A. colBegin(); cur_r != A. rowEnd(); ++ cur_r, ++ cur_c, ++ offset) {
	    
				//try to find the pivot.
				tmp_r = cur_r;

				tmp_cp = cur_c -> begin() + offset;

				while ((tmp_cp != cur_c -> end()) && F.isZero(*tmp_cp)) {
					++ tmp_cp;
					++ tmp_r;
				}

				if (tmp_cp == cur_c -> end()) return 1;

				//if swicth two row if nessary. Each row in dense matrix is stored in contiguous space
				if (tmp_r != cur_r) { 
					P.push_back(std::pair<int,int>(offset, (int)(tmp_cp - cur_c -> begin()) ) );

					std::copy (tmp_r -> begin(), tmp_r -> end(), tmp_v.begin());

					std::copy (cur_r -> begin(), cur_r -> end(), tmp_r -> begin());

					std::copy (tmp_v.begin(), tmp_v.end(), cur_r -> begin());
				}

				// continue gauss elimination
	 
				for(tmp_r = cur_r + 1; tmp_r != A.rowEnd(); ++ tmp_r) {	   
		
					//see if need to update the row
					if (!F.isZero(*(tmp_r -> begin() + offset ))) {
		    
						F.div (tmp_e, *(tmp_r -> begin() + offset), *(cur_r -> begin() + offset));

						F.negin(tmp_e);		    
				
						for ( cur_rp = cur_r ->begin(),tmp_rp =  tmp_r -> begin(); 
						      tmp_rp != tmp_r -> end(); ++ tmp_rp, ++ cur_rp )

							F.axpyin ( *tmp_rp, *cur_rp, tmp_e);

						F.assign(*(tmp_r -> begin() + offset), tmp_e);
					}
				}
	    
			}
	

			//second compute inverse of A.
			DenseMatrix tmp(A);

	
			//2a compute inverse of PA, by solving upper-triangeular system, PA = U^{-1} L.
			typename Field::Element Zero;
			typename Field::Element N_one;
			F.init(Zero,0);
			F.init(N_one, -1);
	
			offset = 0;
			for(cur_c = A.colBegin();cur_c != A. colEnd(); ++ cur_c, ++ offset) {
	    
				for (cur_cp = cur_c -> begin();
				     cur_cp != cur_c -> begin() + offset; ++ cur_cp)
					F.assign (*cur_cp, Zero);

				F.assign(*cur_cp, N_one); ++ cur_cp;

				for (; cur_cp != cur_c -> end(); ++ cur_cp)
					F.negin(*cur_cp);

				//matrix is indexed by 0, instead of 1.	

				for (cur_cp = cur_c -> begin() + (A.rowdim() - 1), tmp_r = tmp.rowBegin() + ( A.rowdim() - 1);
				     cur_cp != cur_c -> begin() - 1; -- cur_cp, -- tmp_r) {
		
					F.assign (tmp_e, *cur_cp);
		
					for(tmp_cp = cur_c -> begin() + (A.rowdim() - 1), tmp_rp = tmp_r -> begin() + ( A.rowdim() -1);
					    tmp_cp != cur_cp; -- tmp_cp, -- tmp_rp) 
						F.axpyin(tmp_e, *tmp_cp, *tmp_rp);
		
					F. div(*cur_cp, tmp_e, *tmp_rp);		

					F.negin(*cur_cp);

				}
	    
			}
	

	
			// 2b, compute inverse of A, A^{-1} = (PA)^{-1} P
			std::vector<std::pair<int, int> >::reverse_iterator v_p;

			for(v_p = P.rbegin(); v_p != P.rend(); ++ v_p) {

				cur_c = A.colBegin() + v_p -> first;
	    
				tmp_c = A.colBegin() + v_p -> second;

				std::copy (cur_c -> begin(), cur_c -> end(), tmp_v.begin());

				std::copy (tmp_c -> begin(), tmp_c -> end(), cur_c -> begin());

				std::copy (tmp_v.begin(), tmp_v.end(), tmp_c -> begin());

			}
         
	
			return 0;
		}
	};

	template<>
	inline long MatrixInverse::matrixInverseIn(const MultiModDouble& F, BlasBlackbox<MultiModDouble>& A) {
		throw LinboxError("LinBox ERROR: use of MultiModDouble with too large moduli is not allowed at this time\n");
		return 0;
	}

} // namespace LinBox


#endif //__LINBOX_matrix_inverse_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
