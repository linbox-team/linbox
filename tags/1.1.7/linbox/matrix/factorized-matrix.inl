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

#ifndef __LINBOX_factorized_matrix_INL
#define __LINBOX_factorized_matrix_INL

namespace LinBox
{


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
				
		FFPACK::applyP( _F, FFLAS::FflasRight, FFLAS::FflasNoTrans, _m,0,_m, L.getWritePointer(), _m, _QQ.getPointer() );
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
	
		FFPACK::applyP( _F, FFLAS::FflasLeft, FFLAS::FflasTrans, _n, 0, _m, S.getWritePointer(), _n, _QQ.getPointer() );
		return S;
	}



	/*
	 * Solvers with Matrices: Operand=BlasMatrix<Element>
	 */

	template <class Field> 
	class FactorizedMatrixLeftSolve<Field, BlasMatrix<typename Field::Element> > {
	public:
		
		BlasMatrix<typename Field::Element>& operator() (const Field& F, 
								 const LQUPMatrix<Field>& A,
								 BlasMatrix<typename Field::Element>& X,
								 const BlasMatrix<typename Field::Element>& B) const{
			linbox_check (A.coldim() == X.rowdim());
			linbox_check (A.rowdim() == B.rowdim());
			linbox_check (B.coldim() == X.coldim());
			int info;
			
			FFPACK::fgetrs (F, FFLAS::FflasLeft, A.rowdim(), A.coldim(), B.coldim(), A.getrank(),
					A.getPointer(), A.getStride(), A.getP().getPointer(), A.getQ().getPointer(),
					X.getPointer(), X.getStride(),
					B.getPointer(), B.getStride(), &info);
			if (info > 0)
				throw LinboxMathInconsistentSystem ("Linear system is inconsistent");
			
			return X;
		}
		
		BlasMatrix<typename Field::Element>& operator() (const Field& F, 
								 const LQUPMatrix<Field>& A,
								 BlasMatrix<typename Field::Element>& B) const{
			
			int info;
			linbox_check (A.coldim() == A.rowdim()); 
			linbox_check (A.coldim() == B.rowdim());
			
			FFPACK::fgetrs (F, FFLAS::FflasLeft, B.rowdim(), B.coldim(), A.getrank(),
					A.getPointer(), A.getStride(),
					A.getP().getPointer(), A.getQ().getPointer(),
					B.getPointer(), B.getStride(), &info);
			if (info > 0)
				throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

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
			linbox_check (A.rowdim() == X.coldim());
			linbox_check (A.coldim() == B.coldim());
			linbox_check (B.rowdim() == X.rowdim());
			int info;
			
			FFPACK::fgetrs (F, FFLAS::FflasRight, A.rowdim(), A.coldim(), B.rowdim(), A.getrank(),
					A.getPointer(), A.getStride(),
					A.getP().getPointer(), A.getQ().getPointer(),
					X.getPointer(), X.getStride(),
					B.getPointer(), B.getStride(), &info);
			if (info > 0)
				throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

			return X;
		}
		
		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const LQUPMatrix<Field>& A,
								  BlasMatrix<typename Field::Element>& B ) const{
			
			int info;
			linbox_check (A.coldim() == A.rowdim()); 
			linbox_check (A.rowdim() == B.coldim());
			
			FFPACK::fgetrs (F, FFLAS::FflasRight, B.rowdim(), B.coldim(), A.getrank(),
					A.getPointer(), A.getStride(),
					A.getP().getPointer(), A.getQ().getPointer(),
					B.getPointer(), B.getStride(), &info);
			if (info > 0)
				throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

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
			linbox_check (A.rowdim() == B.rowdim());
			X = B;
			return  (*this)(F, A, X); 
		}
		
		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const LQUPMatrix<Field>& A,
								  BlasMatrix<typename Field::Element>& B ) const{
		
			linbox_check (A.rowdim() == B.rowdim());

			FFPACK::solveLB2 (F, FFLAS::FflasLeft, B.rowdim(), B.coldim(), A.getrank(),
					  A.getPointer(), A.getStride(), 
					  A.getQ().getPointer(),
					  B.getPointer(), B.getStride());

			return B;
		}
	}; // end of class FactorizedMatrixLeftLSolve
	
	template <class Field> 
	class FactorizedMatrixRightLSolve<Field, BlasMatrix<typename Field::Element> > {
	public:
		BlasMatrix<typename Field::Element>& operator() (const Field& F, 
								 const LQUPMatrix<Field>& A,
								 BlasMatrix<typename Field::Element>& X,
								 const BlasMatrix<typename Field::Element>& B) const{
			linbox_check (A.rowdim() == B.coldim());
			X = B;
			return  (*this)( F, A, X );
		}
		
		BlasMatrix<typename Field::Element>& operator() (const Field& F, 
								 const BlasMatrix<typename Field::element>& A,
								 BlasMatrix<typename Field::Element>& B) const{
			
			linbox_check( A.rowdim() == B.coldim() );

			FFPACK::solveLB2 (F, FFLAS::FflasRight, B.rowdim(), B.coldim(), A.getrank(),
					  A.getPointer(), A.getStride(), 
					  A.getQ().getPointer(), B.getPointer(), B.getStride());
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
			
			linbox_check (A.coldim() == X.rowdim());
			linbox_check (B.coldim() == X.coldim());
			linbox_check (A.rowdim() == B.rowdim());
			
			bool consistent = true;
			size_t ldb = B.getStride();
			size_t ldx = X.getStride();
			typename Field::Element * Bp = B.getPointer();
			typename Field::Element * Xp = X.getPointer();
			typename Field::Element zero,one;
			F.init(zero, 0UL);
			F.init(one, 1UL);
			
			for (size_t i = A.getrank(); i < B.rowdim(); ++i)
				for (size_t j = 0; j < B.coldim(); ++j)
					if (!F.isZero (*(Bp + i*ldb + j)))
						consistent = false;
			if (!consistent) 
				throw  LinboxMathInconsistentSystem ("Linear system is inconsistent");
			
			// The last rows of B are now supposed to be 0

			for (size_t i=0; i < A.getrank(); ++i)
				FFLAS::fcopy (F, B.coldim(), Xp + i*ldx, 1, Bp + i*ldx,1);
			
			FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				      A.getrank(), X.coldim(), one, A.getPointer(), A.getStride(), X.getPointer(), X.getStride());

			for (size_t i=A.getrank(); i < X.rowdim(); ++i)
				for (size_t j = 0; j < X.coldim(); ++j)
					F.assign (*(Xp + i*ldx + j), zero);
			
			return X;
			
		}
		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const BlasMatrix<typename Field::Element>& A, 
								  BlasMatrix<typename Field::Element>& B ) const{
			
			linbox_check (A.coldim() == A.rowdim());
			linbox_check (A.coldim() == B.rowdim());
			typename Field::Element one,zero;
			typename Field::Element * Bp = B.getPointer();
			size_t ldb = B.getStride();
			F.init(one, 1UL);
			F.init(zero, 0UL);
			bool consistent = true;
			
			for (size_t i = A.getrank(); i < B.rowdim(); ++i)
				for (size_t j = 0; j < B.coldim(); ++j)
					if (!F.isZero (*(Bp + i*ldb + j)))
						consistent = false;
			if (!consistent) 
				throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

			FFLAS::ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				      A.getrank(), B.coldim(), one, A.getPointer(), A.getStride(), Bp, ldb);

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
			linbox_check (X.coldim() == A.rowdim());
			linbox_check (X.rowdim() == B.rowdim());
			linbox_check (A.coldim() == B.coldim());
			typename Field::Element one,zero;
			F.init(one, 1UL);
			F.init(zero, 0UL);
			typename Field::Element * Bp = B.getPointer();
			typename Field::Element * Xp = X.getPointer();
			size_t R = A.getrank();
			size_t ldb = B.getStride();
			size_t ldx = X.getStride();
			
			for (size_t i = 0; i < X.getrowdim(); ++i)
				FFLAS::fcopy (F, R, Xp + i*ldx, 1, Bp + i*ldb,1);

			FFLAS::ftrsm (F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
			       X.rowdim(), R, one, A.getPointer(), A.getStride(), X.getPointer(), X.getStride());

			bool consistent = true;
			if (B.coldim() > X.coldim()) {
				typename Field::Element* W = new typename Field::Element [B.rowdim() * (B.coldim() - R)];
				size_t ldw = B.rowdim();
				
				FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
					      X.rowdim(), B.coldim() - R, R,
					      one, Xp, X.getStride(), A.getPointer() + R, A.getStride,
					      zero, W, ldw);
				
				for (size_t i = 0; i < B.rowdim(); ++i)
					for (size_t j = 0; j < B.coldim()-R; ++j)
						if (!F.areEqual (*(W + i*ldw + j), *(Bp + R + i*ldb +j)))
							consistent = false;
				delete[] W;
			} else {
				FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
					      X.rowdim(), B.coldim() - R, R,
					      one, Xp, X.getStride(), A.getPointer() + R, A.getStride,
					      zero, Xp + R, X.getStride());
				
				for (size_t i = 0; i < B.rowdim(); ++i)
					for (size_t j = R; j < B.coldim(); ++j)
						if (!F.areEqual (*(X + i*ldx + j), *(Bp + i*ldb +j)))
							consistent = false;
			}
			if (!consistent)
				throw  LinboxMathInconsistentSystem ("Linear system is inconsistent");

			for (size_t i = 0; i < X.rowdim(); ++i)
				for (size_t j = R; j < X.coldim(); ++j)
					F.assign (*(Xp + i*ldx + j), zero);
			return X;
		}

		BlasMatrix<typename Field::Element>& operator() ( const Field& F, 
								  const LQUPMatrix<Field>& A, 
								  BlasMatrix<typename Field::Element>& B ) const{
			linbox_check (A.coldim() == A.rowdim());
			linbox_check (B.coldim() == A.rowdim());

			typename Field::Element one,zero,mone;
			F.init (one, 1UL);
			F.neg (mone,one);
			F.init (zero, 0UL);
			typename Field::Element * Bp = B.getPointer();
			size_t ldb = B.getStride();
			size_t R = A.getrank();
			
			FFLAS::ftrsm (F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				      B.rowdim(), R, one, A.getPointer(), A.getStride(), B.getPointer(), B.getStride());

			FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				      B.rowdim(), B.coldim() - R, R,
				      one, Bp, B.getStride(), A.getPointer() + R, A.getStride,
				      mone, Bp + R, B.getStride());

			bool consistent = true;
			for (size_t i = 0; i < B.rowdim(); ++i)
				for (size_t j = R; j < B.coldim(); ++j)
					if (!F.isZero (*(Bp + i*ldb + j)))
						consistent = false;
			if (!consistent)
				throw  LinboxMathInconsistentSystem ("Linear system is inconsistent");

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
			linbox_check (A.coldim() == x.size());
			linbox_check (A.rowdim() == b.size());
			int info;

			FFPACK::fgetrs (F, FFLAS::FflasLeft, A.rowdim(), A.coldim(), 1, A.getrank(),
					A.getPointer(), A.getStride(), A.getP().getPointer(), A.getQ().getPointer(),
					&x[0], 1, &b[0], 1, &info);
			if (info > 0)
				throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

			return x;
		}
		
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A, 
								   std::vector<typename Field::Element>& b ) const{
			
			int info;
			linbox_check (A.coldim() == A.rowdim()); 
			linbox_check (A.coldim() == b.size());
			
			FFPACK::fgetrs (F, FFLAS::FflasLeft, b.size(), 1, A.getrank(),
					A.getPointer(), A.getStride(),
					A.getP().getPointer(), A.getQ().getPointer(),
					&b[0], 1, &info);
			if (info > 0)
				throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

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
			linbox_check (A.rowdim() == x.size());
			linbox_check (A.coldim() == b.size());
			int info;
			
			FFPACK::fgetrs (F, FFLAS::FflasRight, A.rowdim(), A.coldim(), 1, A.getrank(),
					A.getPointer(), A.getStride(),
					A.getP().getPointer(), A.getQ().getPointer(),
					&x[0], x.size(), &b[0], b.size(), &info);
			if (info > 0)
				throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

			return x;
		}
		
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A, 
								   std::vector<typename Field::Element>& b ) const{
			
			int info;
			linbox_check (A.coldim() == A.rowdim()); 
			linbox_check (A.rowdim() == b.size());
			
			FFPACK::fgetrs (F, FFLAS::FflasRight, 1, b.size(), A.getrank(),
					A.getPointer(), A.getStride(),
					A.getP().getPointer(), A.getQ().getPointer(),
					&b[0], b.size(), &info);
			if (info > 0)
				throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

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
			linbox_check (A.rowdim() == b.size());

			FFPACK::solveLB2 (F, FFLAS::FflasLeft, b.size(), 1, A.getrank(),
					  A.getPointer(), A.getStride(), 
					  A.getQ().getPointer(), &b[0], 1);

			return b;

			/* BB: unreachable  !
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
			*/
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
			linbox_check (A.rowdim() == b.size());

			FFPACK::solveLB2 (F, FFLAS::FflasRight, 1, b.size(),  A.getrank(),
					  A.getPointer(), A.getStride(), 
					  A.getQ().getPointer(), &b[0], b.size());
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
			linbox_check (A.coldim() == x.size());
			linbox_check (A.rowdim() == b.size());
			
			bool consistent = true;
			typename Field::Element * bp = &b[0];           ;
			typename Field::Element * xp = &x[0];
			typename Field::Element zero,one;
			F.init(zero, 0UL);
			F.init(one, 1UL);
			
			for (size_t i = A.getrank(); i < b.size(); ++i)
				if (!F.isZero (b[i]))
					consistent = false;
			if (!consistent) 
				throw  LinboxMathInconsistentSystem ("Linear system is inconsistent");
			
			// The last rows of B are now supposed to be 0

			for (size_t i=0; i < A.getrank(); ++i)
				F.assign (x[i], b[i]);
			
			FFLAS::ftrsv (F, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				      A.getrank(), A.getPointer(), A.getStride(), xp, 1);

			for (size_t i=A.getrank(); i < x.size(); ++i)
				F.assign (x[i], zero);
			return x;
			
		}
		
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A,
								   std::vector<typename Field::Element>& b ) const{
			
			linbox_check (A.coldim() == A.rowdim());
			linbox_check (A.coldim() == b.size());
			typename Field::Element one,zero;
			F.init(one, 1UL);
			F.init(zero, 0UL);
			bool consistent = true;
			
			for (size_t i = A.getrank(); i < b.size(); ++i)
				if (!F.isZero (b[i]))
					consistent = false;
			if (!consistent) 
				throw  LinboxMathInconsistentSystem ("Linear system is inconsistent");

			FFLAS::ftrsv (F, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 
				      A.getrank(), A.getPointer(), A.getStride(), &b[0], 1);

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
			linbox_check (x.size() == A.rowdim());
			linbox_check (A.coldim() == b.size());
			typename Field::Element one,zero;
			F.init(one, 1UL);
			F.init(zero, 0UL);
			typename Field::Element * bp = b.getPointer();
			typename Field::Element * xp = x.getPointer();
			size_t R = A.getrank();
			
			for (size_t i = 0; i < R; ++i)
				F.assign (x[i], b[i]);
			
			FFLAS::ftrsv (F, FFLAS::FflasUpper, FFLAS::FflasTrans, FFLAS::FflasNonUnit, 
				      R, A.getPointer(), A.getStride(), xp, 1);

			bool consistent = true;
			if (b.size() > x.size()) {
				typename Field::Element* W = new typename Field::Element [b.size() - R];

				FFLAS::fgemv (F, FFLAS::FflasTrans,
					      R, b.size() - R, 
					      one, A.getPointer() + R, A.getStride, xp, 1, 
					      zero, W, 1);
			
				for (size_t i = 0; i < b.size() - R; ++i)
					if (!F.areEqual (W[i], b[i + R]))
						consistent = false;
				delete[] W;
			} else {
				FFLAS::fgemv (F, FFLAS::FflasTrans,
					      R, b.size() - R,
					      one, A.getPointer() + R, A.getStride, xp, 1, 
					      zero, xp + R, 1);
			
				for (size_t i = R; i < b.size(); ++i)
					if (!F.areEqual (x[i], b[i]))
						consistent = false;
			}

			if (!consistent)
				throw  LinboxMathInconsistentSystem ("Linear system is inconsistent");

			for (size_t j = R; j < x.size(); ++j)
				F.assign (x[j], zero);
			return x;
		}
		std::vector<typename Field::Element>& operator() ( const Field& F, 
								   const LQUPMatrix<Field>& A,
								   std::vector<typename Field::Element>& b ) const{
			linbox_check (A.coldim() == A.rowdim());
			linbox_check (b.size() == A.rowdim());

			typename Field::Element one,zero,mone;
			F.init (one, 1UL);
			F.neg (mone,one);
			F.init (zero, 0UL);
			typename Field::Element * bp = &b[0];
			size_t R = A.getrank();
			
			FFLAS::ftrsv (F, FFLAS::FflasUpper, FFLAS::FflasTrans, FFLAS::FflasNonUnit, 
				      R, one, A.getPointer(), A.getStride(), bp, 1);

			FFLAS::fgemv (F, FFLAS::FflasTrans,
				      R, b.size() - R,
				      one, A.getPointer() + R, A.getStride, bp, 1,
				      mone, bp + R, 1);

			bool consistent = true;
			for (size_t j = R; j < b.size(); ++j)
				if (!F.isZero (b[j]))
					consistent = false;
			if (!consistent)
				throw  LinboxMathInconsistentSystem ("Linear system is inconsistent");

			return b;
		}
	}; // end of class FactorizedMatrixRightUSolve


} //end of namespace LinBox


#endif // __LINBOX_factorized_matrix_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
