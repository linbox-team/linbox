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
	template <class Field>
	inline const TriangularBlasMatrix<typename Field::Element>& LQUPMatrix<Field>::getL() const {
		
		TriangularBlasMatrix<typename Field::Element>* L = new TriangularBlasMatrix<typename Field::Element>(_m,_m, low, unit);
		for ( size_t i=0; i<_m; ++i )
			for ( size_t j=0; j<i; ++j )
				L->setEntry( i, j, _L.getEntry(i,j) );
		FFLAPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans, _m,0,_m, L->getPointer(), _m, _Q.getPointer() );
		
		return *L;
		
	}

	// get the matrix U
	template <class Field>
	inline const TriangularBlasMatrix<typename Field::Element>& LQUPMatrix<Field>::getU() const { 

		TriangularBlasMatrix<typename Field::Element>* U = new  TriangularBlasMatrix<typename Field::Element>(_m,_n, up, nonunit);
		for ( size_t i=0; i<_m; ++i )
			for ( size_t j=i; j<_n; ++j )
				U->setEntry( i, j, _U.getEntry(i,j) );
		return *U;
	}

	// get the Matrix S (from the LSP factorization of A deduced from LQUP)
	template <class Field>
	inline const BlasMatrix<typename Field::Element>& LQUPMatrix<Field>::getS() const {
		
		BlasMatrix<typename Field::Element>* S = new BlasMatrix<typename Field::Element>(getU()) ;
		FFLAPACK::applyP( _F, FFLAS::FflasLeft, FFLAS::FflasTrans, _m, 0, _m, S, _m, _Q.getPointer() );
		return *S;
	}



	/*
	 * Solvers with Matrices: Operand=BlasMatrix<Element>
	 */

	template <class Field> 
	class FactorizedMatrixLeftSolve<Field, BlasMatrix<typename Field::Element> > {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::Element>& L, 
				  BlasMatrix<typename Field::Element>& B ) const{
			
			typename Field::Element one;
			F.init( one, 1UL );
			
			// Inversion of L
			// Q = Id since A is invertible
			FFLAS::ftrsm( F, FFLAS::FflasLeft, FFLAS::FflasLower, 
				      FFLAS::FflasNoTrans, FFLAS::FflasUnit, 
				      B.rowdim(), B.coldim(), one,
				      L.getPointer(), L.getStride(), 
				      B.getPointer(), B.getStride() );
			
			// Inversion of U
			FFLAS::ftrsm( F, FFLAS::FflasLeft, FFLAS::FflasUpper, 
				      FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				      B.rowdim(), B.coldim(), one,
				      L.getPointer(), L.getStride(), 
				      B.getPointer(), B.getStride() );
			
			// Inversion of P
			FFLAPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, 
					  B.rowdim(), 0, B.rowdim(), 
					  B.getPointer(), B.getStride(), _P.getPointer() );
			return true;
		}
	}; // end of class FactorizedMatrixLeftSolve

	template <class Field> 
	class FactorizedMatrixRightSolve<Field, BlasMatrix<typename Field::Element> > {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::Element>& L, 
				  BlasMatrix<typename Field::Element>& B ) const{
			
			typename Field::Element one;
			F.init( one, 1UL );
			
			// Inversion of P
			FFLAPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans, 
					  B.coldim(), 0, B.coldim(), B.getPointer(), B.getStride(), _P.getPointer() );
			
			// Inversion of U
			FFLAS::ftrsm( F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				      B.rowdim(), B.coldim(), one,
				      L.getPointer(), L.getStride(), 
				      B.getPointer(), B.getStride() );
			
			// Inversion of L
			// Q = Id since A is invertible
			FFLAS::ftrsm( F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, 
				      B.rowdim(), B.coldim(), one,
				      L.getPointer(), L.getStride(), 
				      B.getPointer(), B.getStride() );
			return true;	
		}
	}; // end of class FactorizedMatrixRightSolve

	template <class Field> 
	class FactorizedMatrixLeftLSolve<Field, BlasMatrix<typename Field::Element> > {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::Element>& L, 
				  BlasMatrix<typename Field::Element>& B, const size_t r ) const{
			size_t m = B.rowdim();
			size_t n = B.coldim();
			if ( m <= n ) {
				FFLAPACK::solveLB( F, FFLAS::FflasLeft, m, n, r, L.getPointer(), L.getStride(), 
						   Q.getPointer(), B.getPointer(), B.getStride() );
			}
			else
				FFLAPACK::solveLB2( F, FFLAS::FflasLeft, m, n, r, L.getPointer(), L.getStride(), 
						    Q.getPointer(), B.getPointer(), B.getStride() );
			return true;
		}
	}; // end of class FactorizedMatrixLeftLSolve
	
	template <class Field> 
	class FactorizedMatrixRightLSolve<Field, BlasMatrix<typename Field::Element> > {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::Element>& L, 
				  BlasMatrix<typename Field::Element>& B, const size_t r ) const{
			size_t m = B.rowdim();
			size_t n = B.coldim();
			if ( m <= n ) {
				FFLAPACK::solveLB( F, FFLAS::FflasRight, m, n, r, L.getPointer(), L.getStride(), 
						   Q.getPointer(), B.getPointer(), B.getStride() );
			}
			else
				FFLAPACK::solveLB2( F, FFLAS::FflasRight, m, n, r, L.getPointer(), L.getStride(), 
						    Q.getPointer(), B.getPointer(), B.getStride() );
			return true;	
		}
	}; // end of class FactorizedMatrixRightLsolve

	template <class Field> 
	class FactorizedMatrixLeftUSolve<Field, BlasMatrix<typename Field::Element> > {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::Element>& U, 
				  BlasMatrix<typename Field::Element>& B ) const{

			FFLAS::ftrsm( F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
				      B.rowdim(), B.coldim(), one, 
				      U.getPointer(), U.getStride(), B.getPointer(), B.getStride() );
		}

	}; // end of class FactorizedMatrixLeftUSolve

	template <class Field> 
	class FactorizedMatrixRightUSolve<Field, BlasMatrix<typename Field::Element> > {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::Element>& U, 
				  BlasMatrix<typename Field::Element>& B ) const{
			FFLAS::ftrsm( F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				      B.rowdim(), B.coldim(), one,
				      U.getPointer(), U.getStride(), B.getPointer(), B.getStride() );	
		}
	}; // end of class FactorizedMatrixRightUSolve


	/*
	 * Solvers with Matrices: Operand=std::vector<Element>
	 */

	template <class Field> 
	class FactorizedMatrixLeftSolve<Field, std::vector<typename Field::Element> > {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::Element>& L, 
				  std::vector<typename Field::Element>& b ) const{
			
			typename Field::Element one;
			F.init( one, 1UL );
			
			// Inversion of L
			// Q = Id since A is invertible
			FFLAS::ftrsv( F, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, 
				      b.size(), L.getPointer(), L.getStride(), &b[0], 1 );
			
			// Inversion of U
			FFLAS::ftrsv( F, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				      b.size(), L.getPointer(), L.getStride(), &b[0], 1 );
			
			// Inversion of P
			FFLAPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, 
					  b.size(), 0, b.size(), 
					  &b[0], b.size(), _P.getPointer() );
			return true;
		}
	}; // end of class FactorizedMatrixLeftSolve

	template <class Field> 
	class FactorizedMatrixRightSolve<Field, std::vector<typename Field::Element> > {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::Element>& L, 
				  std::vector<typename Field::Element>& b ) const{
			
			typename Field::Element one;
			F.init( one, 1UL );
			
			// Inversion of P
			FFLAPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans, 
					  b.size(), 0, b.size(), &b[0], b.size(), _P.getPointer() );
			
			// Inversion of U
			FFLAS::ftrsv( F, FFLAS::FflasUpper, FFLAS::FflasTrans, FFLAS::FflasNonUnit, 
				      b.size(), L.getPointer(), L.getStride(), &b[0], 1 );
			
			// Inversion of L
			// Q = Id since A is invertible
			FFLAS::ftrsv( F, FFLAS::FflasLower, FFLAS::FflasTrans, FFLAS::FflasUnit, 
				      b.size(), L.getPointer(), L.getStride(), &b[0], 1 );
			return true;	
		}
	}; // end of class FactorizedMatrixRightSolve

	template <class Field> 
	class FactorizedMatrixLeftLSolve<Field, std::vector<typename Field::Element> > {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::Element>& L, 
				  std::vector<typename Field::Element>& b, const size_t r ) const{
			size_t m = b.rowdim();
			size_t n = b.coldim();
			if ( m <= n ) {
				FFLAPACK::solveLB( F, FFLAS::FflasLeft, m, n, r, L.getPointer(), L.getStride(), 
						   Q.getPointer(), b.getPointer(), b.getStride() );
			}
			else
				FFLAPACK::solveLB2( F, FFLAS::FflasLeft, m, n, r, L.getPointer(), L.getStride(), 
						    Q.getPointer(), b.getPointer(), b.getStride() );
			return true;
		}
	}; // end of class FactorizedMatrixLeftLSolve
	
	template <class Field> 
	class FactorizedMatrixRightLSolve<Field, std::vector<typename Field::Element> > {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::Element>& L, 
				  std::vector<typename Field::Element>& b, const size_t r ) const{
			size_t m = b.rowdim();
			size_t n = b.coldim();
			if ( m <= n ) {
				FFLAPACK::solveLB( F, FFLAS::FflasRight, m, n, r, L.getPointer(), L.getStride(), 
						   Q.getPointer(), b.getPointer(), b.getStride() );
			}
			else
				FFLAPACK::solveLB2( F, FFLAS::FflasRight, m, n, r, L.getPointer(), L.getStride(), 
						    Q.getPointer(), b.getPointer(), b.getStride() );
			return true;	
		}
	}; // end of class FactorizedMatrixRightLsolve

	template <class Field> 
	class FactorizedMatrixLeftUSolve<Field, std::vector<typename Field::Element> > {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::Element>& U, 
				  std::vector<typename Field::Element>& b ) const{

			FFLAS::ftrsv( F, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
				      bsize(), U.getPointer(), U.getStride(), &b[0], 1 );
		}

	}; // end of class FactorizedMatrixLeftUSolve

	template <class Field> 
	class FactorizedMatrixRightUSolve<Field, std::vector<typename Field::Element> > {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::Element>& U, 
				  std::vector<typename Field::Element>& b ) const{
			FFLAS::ftrsv( F, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				      b.size(), U.getPointer(), U.getStride(), &b[0], 1 );	
		}
	}; // end of class FactorizedMatrixRightUSolve


} //end of namespace LinBox


#endif
