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
	inline TriangularBlasMatrix<typename Field::Element>& LQUPMatrix<Field>::getL(TriangularBlasMatrix<Element>& L) const {
		
		linbox_check( L.coldim() == _m);
		linbox_check( L.rowdim() == _m);
		linbox_check( L.getUpLo() == BlasTag::low);
		linbox_check( L.getDiag() == BlasTag::unit);
	
		typename Field::Element zero,one;
		_F.init( zero, 0UL );
		_F.init( one, 1UL );

		for ( size_t i=0; i<_m; ++i ){
			size_t j=0;
			for (; j< ((i<_n)?i:_n) ; ++j )
				L.setEntry( i, j, _LU.getEntry(i,j) );			
			for (; j<_m; ++j )
				L.setEntry( i, j, zero );
		}		
				
		FFPACK::applyP( _F, FFLAS::FflasRight, FFLAS::FflasNoTrans, _m,0,_m, L.getWritePointer(), _m, _Q.getPointer() );
		for ( size_t i=0; i<_m; ++i )
			L.setEntry( i, i, one );
					
		return L;
		
	}
										      
	// get the matrix U
	template <class Field>
	inline TriangularBlasMatrix<typename Field::Element>& LQUPMatrix<Field>::getU(TriangularBlasMatrix<typename Field::Element>& U) const { 
			
		linbox_check( U.rowdim() == _m);
		linbox_check( U.coldim() == _n);
		linbox_check( U.getUpLo() == BlasTag::up);
		linbox_check( U.getDiag() == BlasTag::nonunit);
		for ( size_t i=0; i<_m; ++i )
			for ( size_t j=i; j<_n; ++j )
				U.setEntry( i, j, _LU.getEntry(i,j) );
		return U;
	}
	
	// get the Matrix S (from the LSP factorization of A deduced from LQUP)
	template <class Field>
	inline BlasMatrix<typename Field::Element>& LQUPMatrix<Field>::getS(BlasMatrix<typename Field::Element>& S) const {
				
		linbox_check( S.rowdim() == _m);
		linbox_check( S.coldim() == _n);  
		typename Field::Element zero;
		_F.init( zero, 0UL);
		for ( size_t i=0; i<_m; ++i ){
			for ( size_t j=0; j<i; ++j )
				S.setEntry( i, j, zero );
			for ( size_t j=i; j<_n; ++j )
				S.setEntry( i, j, _LU.getEntry(i,j) );
		}
	
		FFPACK::applyP( _F, FFLAS::FflasLeft, FFLAS::FflasTrans, _n, 0, _m, S.getWritePointer(), _n, _Q.getPointer() );
		return S;
	}



	/*
	 * Solvers with Matrices: Operand=BlasMatrix<Element>
	 */

	template <class Field> 
	class FactorizedMatrixLeftSolve<Field, BlasMatrix<typename Field::Element> > {
	public:
		
		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const LQUPMatrix<Field>& A,
								  BlasMatrix<typename Field::Element>& X,
								  const BlasMatrix<typename Field::Element>& B ) const{
			linbox_check( A.coldim() == A.rowdim() ); 
			linbox_check( A.coldim() == B.rowdim() );
			linbox_check( A.getrank() == B.rowdim() );
			
			X = B;
			return (*this)( F, A, X );
		}
		
		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const LQUPMatrix<Field>& A,
								  BlasMatrix<typename Field::Element>& B ) const{
			
			size_t m = B.rowdim();
			size_t n = B.coldim();
			linbox_check( A.coldim() == A.rowdim() ); 
			linbox_check( A.coldim() == m );
			linbox_check( A.getrank() == m );
			
			typename Field::Element one, zero;
			F.init( one, 1UL );
			F.init( zero, 0UL);
			
			
			// Inversion of L
			// Q = Id since A is invertible
			FFLAS::ftrsm( F, FFLAS::FflasLeft, FFLAS::FflasLower, 
				      FFLAS::FflasNoTrans, FFLAS::FflasUnit, 
				      m, n, one,
				      A.getPointer(), A.getStride(), 
				      B.getPointer(), B.getStride() );
			
			// Inversion of U
			FFLAS::ftrsm( F, FFLAS::FflasLeft, FFLAS::FflasUpper, 
				      FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				      m, n, one,
				      A.getPointer(), A.getStride(), 
				      B.getPointer(), B.getStride() );
			
			// Inversion of P
			FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, 
					n, 0, m, 
					B.getPointer(), B.getStride(), A.getP().getPointer() );
			
			return B;
		}
	}; // end of class FactorizedMatrixLeftSolve
	
	template <class Field> 
	class FactorizedMatrixRightSolve<Field, BlasMatrix<typename Field::Element> > {
	public:
		
		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const LQUPMatrix<Field>& A,
								  BlasMatrix<typename Field::Element>& X,
								  const BlasMatrix<typename Field::Element>& B ) const{
			
			linbox_check( A.coldim() == A.rowdim() ); 
			linbox_check( A.coldim() == B.rowdim() );
			linbox_check( A.getrank() == B.rowdim() );

			X = B;
			return (*this)( F, A, X );
		}

		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								 const LQUPMatrix<Field>& A, 
								 BlasMatrix<typename Field::Element>& B ) const{
			
			size_t m = B.rowdim();
			size_t n = B.coldim();
			
			linbox_check( A.coldim() == A.rowdim() ); 
			linbox_check( A.rowdim() == n );
			linbox_check( A.getrank() == n );

			typename Field::Element one;
			F.init( one, 1UL );
			
			// Inversion of P
			FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans, 
					  m, 0, n, B.getPointer(), B.getStride(), A.getP().getPointer() );
			
			// Inversion of U
			FFLAS::ftrsm( F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				      m, n, one,
				      A.getPointer(), A.getStride(), 
				      B.getPointer(), B.getStride() );
			
			// Inversion of L
			// Q = Id since A is invertible
			FFLAS::ftrsm( F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, 
				      m, n, one,
				      A.getPointer(), A.getStride(), 
				      B.getPointer(), B.getStride() );
			return B;	
		}
	}; // end of class FactorizedMatrixRightSolve

	template <class Field> 
	class FactorizedMatrixLeftLSolve<Field, BlasMatrix<typename Field::Element> > {
	public:
		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const LQUPMatrix<Field>& A,
								  BlasMatrix<typename Field::Element>& X,
								  const BlasMatrix<typename Field::Element>& B ) const{
			linbox_check( A.rowdim() == B.rowdim() );
			X = B;
			return  (*this)(F, A, X); 
		}
		
		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const LQUPMatrix<Field>& A,
								  BlasMatrix<typename Field::Element>& B ) const{
		
			linbox_check( A.rowdim() == B.rowdim() );
			size_t m = B.rowdim();
			size_t n = B.coldim();
			size_t r = A.getrank();
			if ( A.rowdim() <= A.coldim() ) {
				FFPACK::solveLB( F, FFLAS::FflasLeft, m, n, r, A.getPointer(), A.getStride(), 
						   A.getQ().getPointer(), B.getPointer(), B.getStride() );
			}
			else
				FFPACK::solveLB2( F, FFLAS::FflasLeft, m, n, r, A.getPointer(), A.getStride(), 
						    A.getQ().getPointer(), B.getPointer(), B.getStride() );
			return B;
		}
	}; // end of class FactorizedMatrixLeftLSolve
	
	template <class Field> 
	class FactorizedMatrixRightLSolve<Field, BlasMatrix<typename Field::Element> > {
	public:
		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const LQUPMatrix<Field>& A,
								  BlasMatrix<typename Field::Element>& X,
								  const BlasMatrix<typename Field::Element>& B ) const{
			linbox_check( A.rowdim() == B.coldim() );
			X = B;
			return  (*this)( F, A, X );
		}
		
		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const BlasMatrix<typename Field::element>& A,
								  BlasMatrix<typename Field::Element>& B ) const{
			
			linbox_check( A.rowdim() == B.coldim() );
			size_t m = B.rowdim();
			size_t n = B.coldim();
			size_t r = A.getrank();
			if ( A.rowdim() <= A.coldim() ) {
				FFPACK::solveLB( F, FFLAS::FflasRight, m, n, r, A.getPointer(), A.getStride(), 
						   A.getQ().getPointer(), B.getPointer(), B.getStride() );
			}
			else
				FFPACK::solveLB2( F, FFLAS::FflasRight, m, n, r, A.getPointer(), A.getStride(), 
						    A.getQ().getPointer(), B.getPointer(), B.getStride() );
			return B;	
		}
	}; // end of class FactorizedMatrixRightLsolve

	template <class Field> 
	class FactorizedMatrixLeftUSolve<Field, BlasMatrix<typename Field::Element> > {
	public:
		
		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const LQUPMatrix<Field>& A, 
								  BlasMatrix<typename Field::Element>& X,
								  const BlasMatrix<typename Field::Element>& B ) const{
			
			linbox_check( A.getrank() == B.rowdim() );
			X = B;
			return (*this)( F, A, X );
		}
		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const BlasMatrix<typename Field::Element>& A, 
								  BlasMatrix<typename Field::Element>& B ) const{
			
			linbox_check( A.getrank() == B.rowdim() );
			
			typename Field::Element one;
			F.init( one, 1UL );
			
			FFLAS::ftrsm( F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
				      A.getrank(), B.coldim(), one, 
				      A.getPointer(), A.getStride(), B.getPointer(), B.getStride() );
			return B;
		}
		
	}; // end of class FactorizedMatrixLeftUSolve

	template <class Field> 
	class FactorizedMatrixRightUSolve<Field, BlasMatrix<typename Field::Element> > {
	public:
		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const LQUPMatrix<Field>& A, 
								  BlasMatrix<typename Field::Element>& X, 
								  const BlasMatrix<typename Field::Element>& B ) const{
			linbox_check( A.getrank() == B.coldim() );
			X = B;
			return (*this)( F, A, X );
		}

		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const LQUPMatrix<Field>& A, 
								  BlasMatrix<typename Field::Element>& B ) const{
			linbox_check( A.getrank() == B.coldim() );
			typename Field::Element one;
			F.init( one, 1UL );
			
			FFLAS::ftrsm( F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				      B.rowdim(), A.getrank(), one, A.getPointer(), A.getStride(), B.getPointer(), B.getStride() );	
			return B;
		}
	}; // end of class FactorizedMatrixRightUSolve


	/*
	 * Solvers with vectors: Operand=std::vector<Element>
	 */

	template <class Field> 
	class FactorizedMatrixLeftSolve<Field, std::vector<typename Field::Element> > {
	public:
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A, 
								   std::vector<typename Field::Element>& x,
								   const std::vector<typename Field::Element>& b ) const{
			//linbox_check( A.coldim() == A.rowdim() );
			//linbox_check( A.rowdim() == b.size() );
			//linbox_check( A.getrank() == A.rowdim() );
			
			x = b;
			return (*this)( F, A, x );
		}
		
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A, 
								   std::vector<typename Field::Element>& b ) const{
			
			
			size_t m, n, r;
			m = A.rowdim();
			n = A.coldim();
			r = A.getrank(); 
			typename Field::Element zero;
			F.init(zero, 0UL);
			linbox_check( m <= n ); 
							       
			typename Field::Element one;
			F.init( one, 1UL );
			
			if ((r==m)&&(m==n)) {
				linbox_check( m == b.size() );
				
				// Inversion of L
				// Q = Id since A is invertible
				FFLAS::ftrsv( F, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, 
					      b.size(), A.getPointer(), A.getStride(), &b[0], 1 );
				
				// Inversion of U
				FFLAS::ftrsv( F, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
					      b.size(), A.getPointer(), A.getStride(), &b[0], 1 );
							
				// Inversion of P
				FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, 
						1, 0, b.size(), &b[0], 1, A.getP().getPointer() );				
			}
			else {	
				b.resize(n,zero);
		
				// b:= L^(-1).b 
				FFLAS::ftrsv( F, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, 
					      m, A.getPointer(), A.getStride(), &b[0], 1 );
				
				bool consistent=true;					
				if (r < m){
					// b:= Q^t.b
					FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, 
							1, 0, m, &b[0], 1, A.getQ().getPointer() );
					
					// check consistency
					for (size_t i=r;i<m;++i)
						if (!F.isZero(b[i])){
							consistent=false;
							break;
						}
				}
			
				if (consistent){
					// b:= Ur^(-1).b with Ur the rxr leading minor of U
					FFLAS::ftrsv( F, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
						      r, A.getPointer(), A.getStride(), &b[0], 1 );
														
					// b:= P^t.b
					FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, 
							1, 0, n, &b[0], 1, A.getP().getPointer() );
				}
				else{
					for (size_t i=0;i<n;++i)
						F.assign(b[i], zero);	
				}				
			}
			
			return b;
		}
	}; // end of class FactorizedMatrixLeftSolve

	template <class Field> 
	class FactorizedMatrixRightSolve<Field, std::vector<typename Field::Element> > {
	public:
			
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A, 
								   std::vector<typename Field::Element>& x,
								   const std::vector<typename Field::Element>& b ) const{
			linbox_check( A.coldim() == A.rowdim() ); 
			linbox_check( A.rowdim() == b.size() );
			linbox_check( A.getrank() == A.rowdim() );
			
			x = b;
			return (*this)( F, A, x );
		}
		
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A, 
								   std::vector<typename Field::Element>& b ) const{
			
			typename Field::Element one;
			F.init( one, 1UL );
			
			// Inversion of P
			FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans, 
					  1, 0, b.size(), &b[0], 1, A.getP().getPointer() );
			
			// Inversion of U
			FFLAS::ftrsv( F, FFLAS::FflasUpper, FFLAS::FflasTrans, FFLAS::FflasNonUnit, 
				      b.size(), A.getPointer(), A.getStride(), &b[0], 1 );
			
			// Inversion of L
			// Q = Id since A is invertible
			FFLAS::ftrsv( F, FFLAS::FflasLower, FFLAS::FflasTrans, FFLAS::FflasUnit, 
				      b.size(), A.getPointer(), A.getStride(), &b[0], 1 );
			return b;	
		}
	}; // end of class FactorizedMatrixRightSolve

	template <class Field> 
	class FactorizedMatrixLeftLSolve<Field, std::vector<typename Field::Element> > {
	public:
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A, 
								   std::vector<typename Field::Element>& x,
								   const std::vector<typename Field::Element>& b ) const{
			linbox_check( A.rowdim() == b.size() );
			x = b;
			return (*this)( F, A, x );
		}
		
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A,
								   std::vector<typename Field::Element>& b ) const{
			size_t n = b.size(); // bds: b not B
			linbox_check( A.rowdim() == n );
			size_t r = A.getrank();
			
			// To be changed: solveLB is designed for matrices, not for vectors
			if ( A.rowdim() <= A.coldim() ) {
				FFPACK::solveLB( F, FFLAS::FflasLeft, n, 1, r, A.getPointer(), A.getStride(), 
						   A.getQ().getPointer(), &b[0], b.size() );
			}
			else
				FFPACK::solveLB2( F, FFLAS::FflasLeft, n, 1, r, A.getPointer(), A.getStride(), 
						    A.getQ().getPointer(), b.getPointer(), b.getStride() );
			return b;
		}
	}; // end of class FactorizedMatrixLeftLSolve
	
	template <class Field> 
	class FactorizedMatrixRightLSolve<Field, std::vector<typename Field::Element> > {
	public:
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A, 
								   std::vector<typename Field::Element>& x,
								   const std::vector<typename Field::Element>& b ) const{
			linbox_check( A.rowdim() == b.size() );
			x = b;
			return (*this)( F, A, x );
		}
		
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A,
								   std::vector<typename Field::Element>& b ) const{
			size_t n = b.size();
			linbox_check( A.rowdim() == n );
			size_t r = A.getrank();
			
			// To be changed: solveLB is designed for matrices, not for vectors
			if ( A.rowdim() <= A.coldim() ) {
				FFPACK::solveLB( F, FFLAS::FflasRight, 1, n, r, A.getPointer(), A.getStride(), 
						   A.getQ().getPointer(), b.getPointer(), b.getStride() );
			}
			else
				FFPACK::solveLB2( F, FFLAS::FflasRight, 1, n, r, A.getPointer(), A.getStride(), 
							//bds: A.getQ not getQ
						    A.getQ().getPointer(), b.getPointer(), b.getStride() );
			return b;	
		}
	}; // end of class FactorizedMatrixRightLsolve

	template <class Field> 
	class FactorizedMatrixLeftUSolve<Field, std::vector<typename Field::Element> > {
	public:
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A, 
								   std::vector<typename Field::Element>& x,
								   const std::vector<typename Field::Element>& b ) const{
			linbox_check( A.getrank() == b.size() );
			x = b;
			return (*this)( F, A, x );
		}
		
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A,
								   std::vector<typename Field::Element>& b ) const{
			
			linbox_check( A.getrank() == b.size() );
			FFLAS::ftrsv( F, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
				      b.size(), A.getPointer(), A.getStride(), &b[0], 1 );
			//return b, by Z. W
			return b;
		}

	}; // end of class FactorizedMatrixLeftUSolve

	template <class Field> 
	class FactorizedMatrixRightUSolve<Field, std::vector<typename Field::Element> > {
	public:
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A, 
								   std::vector<typename Field::Element>& x,
								   const std::vector<typename Field::Element>& b ) const{
			linbox_check( A.getrank() == b.size() );
			x = b;
			return (*this)( F, A, x );
		}
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A,
								   std::vector<typename Field::Element>& b ) const{
			linbox_check( A.rowdim() == b.size() );
			FFLAS::ftrsv( F, FFLAS::FflasUpper, FFLAS::FflasTrans, FFLAS::FflasNonUnit, 
				      b.size(), A.getPointer(), A.getStride(), &b[0], 1 );
			return b;
		}
	}; // end of class FactorizedMatrixRightUSolve


} //end of namespace LinBox


#endif
