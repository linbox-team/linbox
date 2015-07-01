/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/lsp-tools.h
 * Copyright (C) 2003 Pascal Giorgi
 *
 * Written by Pascal Giorgi pascal.giorgi@ens-lyon.fr
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

#ifndef __LSP_TOOLS_H
#define __LSP_TOOLS_H

#include <vector>
#include "linbox/FFLAS/fflas.h"


/* 
 * Definition of useful function involved in the LSP decomposition
 */


namespace LinBox{


	/* Application of a column permutation.
	 *
	 * F is the field of the computation.
	 * A is the pointer representing the matrix.
	 * m is the row size of the matrix.
	 * n is the column size of the matrix.
	 * lda is the stride of the storage of Elements.
	 * P is a vector coding the permutation.
	 */

	template <class Field>
	void ApplyColPerm (const Field& F,
			   typename Field::Element* A, int m, int n , int lda ,
			   const std::vector<int>& P) {
	
		typename Field::Element copy[n];
		for (int i=0;i<m;i++) {
			for (int j=0;j<n;j++)
				F.assign( copy[j] , *(A+i*lda+j));  
		
			for (int k=0;k<n;k++)
				F.assign( *(A+i*lda+P[k]) , copy[k]);
		}
	
	}
	

	/* Application of a row permutation.
	 *
	 * F is the field of the computation.
	 * A is the pointer representing the matrix.
	 * m is the row size of the matrix.
	 * n is the column size of the matrix.
	 * lda is the stride of the storage of Elements.
	 * P is a vector coding the permutation.
	 */

	template <class Field>
	void ApplyRowPerm (const Field& F,
			   typename Field::Element* A, int m, int n , int lda ,
			   const std::vector<int>& P) {
	
		typename Field::Element copy[m];
		for (int j=0;j<n;j++){
			for (int i=0;i<m;i++)
				F.assign( copy[i] , *(A+i*lda+j));  
		
			for (int k=0;k<m;k++)
				F.assign( *(A+j+P[k]*lda) , copy[k]);
		}
	
	}


	/* Application of a transposed of column permutation.
	 *
	 * F is the field of the computation.
	 * A is the pointer representing the matrix.
	 * m is the row size of the matrix.
	 * n is the column size of the matrix.
	 * lda is the stride of the storage of Elements.
	 * P is a vector coding the permutation.
	 */
	template <class Field>
	void ApplyColPermTrans (const Field& F,
				typename Field::Element* A, int m, int n , int lda ,
				const std::vector<int>& P) {
	
		typename Field::Element copy[n];
		for (int i=0;i<m;i++) {
			for (int j=0;j<n;j++)
				F.assign( copy[j] , *(A+i*lda+j));  

			for (int k=0;k<n;k++)			
				F.assign( *(A+i*lda+k) , copy[P[k]]);
		}
	
	}

	/* Application of a transposed of row permutation.
	 *
	 * F is the field of the computation.
	 * A is the pointer representing the matrix.
	 * m is the row size of the matrix.
	 * n is the column size of the matrix.
	 * lda is the stride of the storage of Elements.
	 * P is a vector coding the permutation.
	 */

	template <class Field>
	void ApplyRowPermTrans (const Field& F,
				typename Field::Element* A, int m, int n , int lda ,
				const std::vector<int>& P) {
		
		typename Field::Element copy[m];
		for (int j=0;j<n;j++){
			for (int i=0;i<m;i++)
				F.assign( copy[i] , *(A+i*lda+j));  
			
			for (int k=0;k<m;k++)
				F.assign( *(A+j+k*lda) , copy[P[k]]);
		}
		
	}


	// this function compute G as it is described in (Bini & Pan - Polynomial and Matrix computation - LSP FACTORS  p.103)
	template <class Field>
	void ComputeG (const Field& F,
		       typename Field::Element* S, int m, int r, int lds,
		       typename Field::Element* A, int mA, int lda,
		       typename Field::Element* L, int ldl) {
	
		typedef typename Field::Element Element;

		std::vector<int> perm(m);
		Element T[r*r];	
		int idx=0;
		int idz=r;
		for (int i=0;i<m;i++)
			if ( !(F.isZero(*(S+idx+i*lds)))) {
				for (int k=idx;k<r;k++)
					F.assign( *(T+idx*r+k), *(S+k+i*lds));
				perm[i]=idx;
				idx++;
			}
			else {
				perm[i]=idz;
				idz++;
			}
		
		// L= A*T^-1
		Field_trsm_up_right (F,mA,r,A,lda,T,r,L,ldl);
		
		// L= L * Perm.
		if (m != r)
			ApplyColPermTrans (F,L,mA,m,ldl,perm);
		
	}
		
	
				
		   
} // end of namespace LinBox 

#endif
