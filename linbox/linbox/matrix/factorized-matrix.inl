/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/matrix/factorized-matrix.inl
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
 *               Clément Pernet clement.pernet@imag.fr
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

#ifndef __FACTORIZED_MATRIX_INL
#define __FACTORIZED_MATRIX_INL

namespace LinBox{


	// get the Matrix L
	template <class Field,class Matrix>
	inline const BlasMatrix<Matrix>& LQUPMatrix<Field,Matrix>::getL() const{


	}

	// get the matrix U
	template <class Field,class Matrix>
	inline const BlasMatrix<Matrix>& LQUPMatrix<Field,Matrix>::getU() const { 

		Element zero;
		F.init(zero,0L);
		BlasMatrix<Matrix> tmp(_LU);
		for (int i=0;i<_m;++i)
			for (int j=0;j<i;++j)
				tmp.setEntry(i,j,zero);			
		return *(new BlasMatrix<Matrix> (tmp));
	}



	/*
	 * Solvers with Matrices
	 */
	// solve AX=B
	template <class Field,class Matrix>
	inline bool LQUPMatrix<Field,Matrix>::solve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const{}

	// solve AX=B (X is stored in B)
	template <class Field,class Matrix>
	inline bool LQUPMatrix<Field,Matrix>::solve(BlasMatrix<Matrix>& B) const{}

	// solve LX=B (L from LQUP)
	template <class Field,class Matrix>
	inline bool LQUPMatrix<Field,Matrix>::Lsolve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const{}
		
	// solve LX=B (L from LQUP) (X is stored in B)
	template <class Field,class Matrix>
	inline bool LQUPMatrix<Field,Matrix>::Lsolve(BlasMatrix<Matrix>& B) const{}

	// solve XU=B (U from LQUP)
	template <class Field,class Matrix>
	inline bool LQUPMatrix<Field,Matrix>:: Usolve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const{}
		
	// solve XU=B (U from LQUP) (X is stored in B)
	template <class Field,class Matrix>
	inline bool LQUPMatrix<Field,Matrix>::Usolve(BlasMatrix<Matrix>& B) const{}


	/*
	 * Solvers with vectors
	 */
	// solve Ax=b
	template <class Field,class Matrix>
	inline bool LQUPMatrix<Field,Matrix>::solve(std::vector<Element>& x, const std::vector<Element>& b) const{}
		
	// solve Ax=b (x is stored in b)
	template <class Field,class Matrix>
	inline bool LQUPMatrix<Field,Matrix>::solve(std::vector<Element>& b) const{}
	
	// solve Lx=b (L from LQUP) 
	template <class Field,class Matrix>
	inline bool LQUPMatrix<Field,Matrix>::Lsolve(std::vector<Element>& x, const std::vector<Element>& b) const{}
				
	// solve Lx=b (L from LQUP) (x is stored in b)
	template <class Field,class Matrix>
	inline bool LQUPMatrix<Field,Matrix>::Lsolve(std::vector<Element>& b) const{}		

	// solve xU=b (U from LQUP)
	template <class Field,class Matrix>
	inline bool LQUPMatrix<Field,Matrix>::Usolve(std::vector<Element>& x, const std::vector<Element>& b) const{}
				
	// solve xU=b (U from LQUP) (x is tored in b)
	template <class Field,class Matrix>
	inline bool LQUPMatrix<Field,Matrix>::Usolve(std::vector<Element>& b) const{}



} //end of namespace LinBox;


#endif
