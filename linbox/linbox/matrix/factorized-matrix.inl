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
	inline const TriangularBlasMatrix<Matrix>& LQUPMatrix<Field,Matrix>::getL() const {
		
		TriangularBlasMatrix<Matrix>* L = new  TriangularBlasMatrix<Matrix>(_m,_m, low, unit);
		for ( size_t i=0; i<_m; ++i )
			for ( size_t j=0; j<i; ++j )
				L->setEntry( i, j, _L.getEntry(i,j) );
		FFLAPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans, _m,0,_m, L->getPointer(), _m, _Q.getPointer() );
		
		return *L;
		
	}

	// get the matrix U
	template <class Field,class Matrix>
	inline const TriangularBlasMatrix<Matrix>& LQUPMatrix<Field,Matrix>::getU() const { 

		TriangularBlasMatrix<Matrix>* U = new  TriangularBlasMatrix<Matrix>(_m,_n, up, nonunit);
		for ( size_t i=0; i<_m; ++i )
			for ( size_t j=i; j<_n; ++j )
				U->setEntry( i, j, _U.getEntry(i,j) );
		return *U;
	}

	// get the Matrix S (from the LSP factorization of A deduced from LQUP)
	template <class Field,class Matrix>
	inline const BlasMatrix<Matrix>& LQUPMatrix<Field,Matrix>::getS() const {
		
		BlasMatrix<Matrix> S( getU() );
		FFLAPACK::applyP( _F, FFLAS::FflasLeft, FFLAS::FflasTrans, _m, 0, _m, S, _m, _Q.getPointer() );
	}



	/*
	 * Solvers with Matrices
	 */
	// solve AX=B
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::left_solve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const{
		
		if ( ( _m != _n ) || ( _rank != _n ) || ( B.coldim() != _n) )
			return false;
		X = B;
		return left_solve ( X );

	}
	
	// solve AX=B (X is stored in B)
	template<class Field> 
	inline bool LQUPMatrix::left_solve(BlasMatrix<typename Field::Element>& B) const{
		
		typename Field::Element one;
		_F.init( one, 1UL );
		
		if ( ( _m != _n ) || ( _rank != _n ) || ( B.coldim() != _n) )
			return false;
		// Inversion of L
		// Q = Id since A is invertible
		FFLAS::ftrsm( _F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, 
			      B.rowdim(), B.coldim(), one,
			      _LU.getPointer(), _LU.getStride(), 
			      B.getPointer(), B.getStride() );
		
		// Inversion of U
		FFLAS::ftrsm( _F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
			      B.rowdim(), B.coldim(), one,
			      _LU.getPointer(), _LU.getStride(), 
			      B.getPointer(), B.getStride() );

		// Inversion of P
		FFLAPACK::applyP( _F, FFLAS::FflasLeft, FFLAS::FflasTrans, _m, 0, _m, B.getPointer(), _m, _P.getPoiner() );
		
	}

	// solve XA=B
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::right_solve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const{	

		if ( ( _m != _n ) || ( _rank != _n ) || ( B.rowldim() != _m) )
			return false;
		X = B;
		return right_solve ( X );
	}

	// solve XA=B (X is stored in B)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::right_solve(BlasMatrix<Matrix>& B) const{
		typename Field::Element one;
		_F.init( one, 1UL );
		
		if ( ( _m != _n ) || ( _rank != _n ) || ( B.coldim() != _n) )
			return false;
		
		// Inversion of P
		FFLAPACK::applyP( _F, FFLAS::FflasLeft, FFLAS::FflasTrans, _m, 0, _m, B.getPointer(), _m, _P.getPoiner() );

		// Inversion of U
		FFLAS::ftrsm( _F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
			      B.rowdim(), B.coldim(), one,
			      _LU.getPointer(), _LU.getStride(), 
			      B.getPointer(), B.getStride() );
		
		// Inversion of L
		// Q = Id since A is invertible
		FFLAS::ftrsm( _F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, 
			      B.rowdim(), B.coldim(), one,
			      _LU.getPointer(), _LU.getStride(), 
			      B.getPointer(), B.getStride() );
	}


	// solve LX=B (L from LQUP)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::left_Lsolve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const{}
		
	// solve LX=B (L from LQUP) (X is stored in B)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::left_Lsolve(BlasMatrix<Matrix>& B) const{}

	// solve XL=B (L from LQUP)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::right_Lsolve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const{}
		
	// solve XL=B (L from LQUP) (X is stored in B)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::right_Lsolve(BlasMatrix<Matrix>& B) const{}


	// solve UX=B (U from LQUP)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::left_Usolve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const{}
		
	// solve UX=B (U from LQUP) (X is stored in B)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::rleft_Usolve(BlasMatrix<Matrix>& B) const{}

	// solve XU=B (U from LQUP)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::right_Usolve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const{}
		
	// solve XU=B (U from LQUP) (X is stored in B)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::right_Usolve(BlasMatrix<Matrix>& B) const{}


	/*
	 * Solvers with vectors
	 */
		
	// solve Ax=b
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::left_solve(std::vector<Element>& x, const std::vector<Element>& b) const{}
		
	// solve Ax=b (x is stored in b)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::left_solve(std::vector<Element>& b) const{}

	// solve xA=b
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::right_solve(std::vector<Element>& x, const std::vector<Element>& b) const{}

	// solve xA=b (x is stored in b)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::right_solve(std::vector<Element>& b) const{}


	// solve Lx=b (L from LQUP)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::left_Lsolve(std::vector<Element>& x, const std::vector<Element>& b) const{}
		
	// solve Lx=b (L from LQUP) (x is stored in b)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::left_Lsolve(std::vector<Element>& b) const{}

	// solve xL=b (L from LQUP)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::right_Lsolve(std::vector<Element>& x, const std::vector<Element>& b) const{}
		
	// solve xL=b (L from LQUP) (x is stored in b)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::right_Lsolve(std::vector<Element>& b) const{}


	// solve Ux=b (U from LQUP)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::left_Usolve(std::vector<Element>& x, const std::vector<Element>& b) const{}
		
	// solve Ux=b (U from LQUP) (x is stored in b)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::rleft_Usolve(std::vector<Element>& b) const{}

	// solve xU=b (U from LQUP)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::right_Usolve(std::vector<Element>& x, const std::vector<Element>& b) const{}
		
	// solve xU=b (U from LQUP) (x is stored in b)
	template<class Field, class Matrix> 
	inline bool LQUPMatrix::right_Usolve(std::vector<Element>& b) const{}



} //end of namespace LinBox


#endif
