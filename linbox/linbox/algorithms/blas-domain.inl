/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/blas-domain.inl
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



#ifndef __BLAS_MATRIX_DOMAIN_INL
#define __BLAS_MATRIX_DOMAIN_INL

#include <linbox/matrix/blas-matrix.h>
#include <linbox/matrix/factorized-matrix.h>

namespace LinBox {


	/*
	 * **********************************************
	 * *** Specialization for BlasMatrix<Element> ***
	 * **********************************************
	 */	


	// Inversion
	template <class Field>	
	class BlasMatrixDomainInv<Field,BlasMatrix<typename Field::Element> > {
	public:
		BlasMatrix<typename Field::Element>& operator() (const Field& F, 
								 BlasMatrix<typename Field::Element>& Ainv,
								 const BlasMatrix<typename Field::Element>& A) const{
			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.rowdim() == Ainv.rowdim());
			linbox_check( A.coldim() == Ainv.coldim());
			
			BlasMatrix<typename Field::Element> tmp(A);			
			return (*this)(F,Ainv,tmp);
		}

		BlasMatrix<typename Field::Element>& operator() (const Field& F, 
								 BlasMatrix<typename Field::Element>& Ainv,
								 BlasMatrix<typename Field::Element>& A) const{
			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.rowdim() == Ainv.rowdim());
			linbox_check( A.coldim() == Ainv.coldim());
			
			FFLAPACK::Invert(F,A.rowdim(),A.getPointer(),A.getStride(),Ainv.getPointer(),Ainv.getStride());
			return Ainv;
		}
		
	};

	// Rank
	template <class Field>
	class 	BlasMatrixDomainRank<Field,BlasMatrix<typename Field::Element> > {
	public:
		inline unsigned int operator() (const Field& F,const BlasMatrix<typename Field::Element>& A) const{
			BlasMatrix<typename Field::Element> tmp(A);
			return (*this)(F,tmp);
		}	

		inline unsigned int 
		operator() (const Field& F, BlasMatrix<typename Field::Element>& A) const{
			return FFLAPACK::Rank(F, A.rowdim(), A.coldim(),A.getPointer(), A.getStride());
		}
	};

	// determinant
	template <class Field>
	class BlasMatrixDomainDet<Field,BlasMatrix<typename Field::Element> > {
	public:
		inline typename Field::Element operator()(const Field& F,const BlasMatrix<typename Field::Element>& A) const{
			BlasMatrix<typename Field::Element> tmp(A);
			return  (*this)(F,tmp);
		}

		inline typename Field::Element operator() (const Field& F,BlasMatrix<typename Field::Element>& A) const{
			return FFLAPACK::Det(F, A.rowdim(), A.coldim(),A.getPointer(), A.getStride());
		}
	};


	/*
	 * specialization for Operand1, Operand2 and Operand3  of type BlasMatrix<Element>
	 */
	
	//  general matrix-matrix multiplication and addition with scaling
	// D= beta.C + alpha.A*B
	template<class Field>
	class 	BlasMatrixDomainMulAdd<Field,BlasMatrix<typename Field::Element>,BlasMatrix<typename Field::Element>, BlasMatrix<typename Field::Element> > {
	public:
		BlasMatrix<typename Field::Element>& operator()(const Field& F,
								BlasMatrix<typename Field::Element>& D, 
								const typename Field::Element& beta, 
								const BlasMatrix<typename Field::Element>& C,
								const typename Field::Element& alpha, 
								const BlasMatrix<typename Field::Element>& A, 
								const BlasMatrix<typename Field::Element>& B) const{
			linbox_check( A.coldim() == B.rowdim());
			linbox_check( C.rowdim() == A.rowdim());
			linbox_check( C.coldim() == B.coldim());
			linbox_check( D.rowdim() == C.rowdim());
			linbox_check( D.coldim() == C.coldim());
			
			D=C;
			
			FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				      C.rowdim(), C.coldim(), A.coldim(),
				      alpha,
				      A.getPointer(), A.getStride(),
				      B.getPointer(), B.getStride(),
				      beta,
				      D.getPointer(), D.getStride());
			return D;
		}

		
		BlasMatrix<typename Field::Element>&operator() (const Field& F,
								const typename Field::Element& beta,
								BlasMatrix<typename Field::Element>& C,
								const typename Field::Element& alpha, 
								const BlasMatrix<typename Field::Element>& A, 
								const BlasMatrix<typename Field::Element>& B) const{
			linbox_check( A.coldim() == B.rowdim());
			linbox_check( C.rowdim() == A.rowdim());
			linbox_check( C.coldim() == B.coldim());
			
			FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				      C.rowdim(), C.coldim(), A.coldim(),
				      alpha,
				      A.getPointer(), A.getStride(),
				      B.getPointer(), B.getStride(),
				      beta,
				      C.getPointer(), C.getStride());
			return C;
		}
	};



	/*
	 * specialization for Operand1 and Operand3 of type std::vector<Element>
	 * and Operand2 of type BlasMatrix<Element>
	 */

	//  general matrix-vector multiplication and addition with scaling
	// d = beta.c + alpha.A*b
	template<class Field>
	class BlasMatrixDomainMulAdd<Field,std::vector<typename Field::Element>,BlasMatrix<typename Field::Element>,std::vector<typename Field::Element> > {
	public:
		std::vector<typename Field::Element>& operator() (const Field& F,
								  std::vector<typename Field::Element>& d, 
								  const typename Field::Element& beta, 
								  const std::vector<typename Field::Element>& c,
								  const typename Field::Element& alpha, 
								  const BlasMatrix<typename Field::Element>& A, 
								  const std::vector<typename Field::Element>& b) const{
			linbox_check( A.coldim() == b.size());
			linbox_check( c.size()   == b.size());
			linbox_check( d.size()   == c.size());
			d=c;
			
			FFLAS::fgemv( F, FFLAS::FflasNoTrans, 
				      A.rowdim(), A.coldim(),
				      alpha,
				      A.getPointer(), A.getStride(),
				      &b[0],1,
				      beta,
				      &d[0],1);
			return d;
		}
		

		std::vector<typename Field::Element>& operator() (const Field& F,
								  const typename Field::Element& beta,
								  std::vector<typename Field::Element>& c,
								  const typename Field::Element& alpha, 
								  const BlasMatrix<typename Field::Element>& A, 
								  const std::vector<typename Field::Element>& b) const{
			linbox_check( A.coldim() == b.size());
			linbox_check( c.size()   == b.size());
			
			FFLAS::fgemv( F, FFLAS::FflasNoTrans, 
				      A.rowdim(), A.coldim(),
				      alpha,
				      A.getPointer(), A.getStride(),
				      &b[0],1,
				      beta,
				      &c[0],1);
			return c;
		}
	};

	template<class Field>
	class BlasMatrixDomainMulAdd<Field,std::vector<typename Field::Element>,std::vector<typename Field::Element>,BlasMatrix<typename Field::Element> > {
	public:
		std::vector<typename Field::Element>& operator() (const Field& F,
								  std::vector<typename Field::Element>& d, 
								  const typename Field::Element& beta, 
								  const std::vector<typename Field::Element>& c,
								  const typename Field::Element& alpha, 
								  const std::vector<typename Field::Element>& a, 
								  const BlasMatrix<typename Field::Element>& B) const{
			linbox_check( B.rowdim() == a.size());
			linbox_check( c.size()   == a.size());
			linbox_check( d.size()   == c.size());
			d=c;
			
			FFLAS::fgemv( F, FFLAS::FflasTrans, 
				      B.rowdim(), B.coldim(),
				      alpha,
				      B.getPointer(), B.getStride(),
				      &a[0],1,
				      beta,
				      &d[0],1);
			return d;
		}
		

		std::vector<typename Field::Element>& operator() (const Field& F,
								  const typename Field::Element& beta,
								  std::vector<typename Field::Element>& c,
								  const typename Field::Element& alpha, 
								  const std::vector<typename Field::Element>& a, 
								  const BlasMatrix<typename Field::Element>& B) const{
			linbox_check( B.rowdim() == a.size());
			linbox_check( c.size()   == a.size());
			
			FFLAS::fgemv( F, FFLAS::FflasTrans, 
				      B.rowdim(), B.coldim(),
				      alpha,
				      B.getPointer(), B.getStride(),
				      &a[0],1,
				      beta,

				      &c[0],1);
			return c;
		}
	};


	/*
	 * Specialization for Operand of type BlasMatrix<Element>
	 */

	template <class Field>
	class BlasMatrixDomainLeftSolve<Field,BlasMatrix<typename Field::Element>,BlasMatrix<typename Field::Element> > {
	public:
		BlasMatrix<typename Field::Element>& operator() (const Field& F,
								 BlasMatrix<typename Field::Element>& X,
								 const BlasMatrix<typename Field::Element>& A,
								 const BlasMatrix<typename Field::Element>& B) const {
			LQUPMatrix<Field> LQUP(F,A);
			LQUP.left_solve(X,B);
			return X;
		}
	

		BlasMatrix<typename Field::Element>& operator() (const Field& F,
								 const BlasMatrix<typename Field::Element>& A, 
								 BlasMatrix<typename Field::Element>& B) const{
			LQUPMatrix<Field> LQUP(F,A);
			LQUP.left_solve(B);
			return B;
		}
	};

	template <class Field>
	class BlasMatrixDomainRightSolve<Field,BlasMatrix<typename Field::Element>,BlasMatrix<typename Field::Element> > {
	public:
		BlasMatrix<typename Field::Element>& operator() (const Field& F,
								 BlasMatrix<typename Field::Element>& X,
								 const BlasMatrix<typename Field::Element>& A,
								 const BlasMatrix<typename Field::Element>& B) const{
			LQUPMatrix<Field> LQUP(F,A);
			LQUP.right_solve(X,B);
			return X;
		}
	
	
		BlasMatrix<typename Field::Element>& operator() (const Field& F,
								 const BlasMatrix<typename Field::Element>& A, 
								 BlasMatrix<typename Field::Element>& B) const{
			LQUPMatrix<Field> LQUP(F,A);
			LQUP.right_solve(B);
			return B;
		}
	
	};

	/*
	 * Specialization for Operand of type std::vector<Element>
	 */

	template <class Field>
	class BlasMatrixDomainLeftSolve<Field, std::vector<typename Field::Element>, BlasMatrix<typename Field::Element> > {
	public:
		std::vector<typename Field::Element>& operator() (const Field& F,
								  std::vector<typename Field::Element>& X,
								  const BlasMatrix<typename Field::Element>& A,
								  const std::vector<typename Field::Element>& B) const {
			LQUPMatrix<Field> LQUP(F,A);
			LQUP.left_solve(X,B);
			return X;
		}
	
		std::vector<typename Field::Element>& operator()(const Field& F,
								 const BlasMatrix<typename Field::Element>& A, 
								 std::vector<typename Field::Element>& B) const{
			LQUPMatrix<Field> LQUP(F,A);
			LQUP.left_solve(B);
			return B;	
		}
		
	};
	
	template <class Field>
	class BlasMatrixDomainRightSolve<Field, std::vector<typename Field::Element>, BlasMatrix<typename Field::Element> > {
	public:
		std::vector<typename Field::Element>& operator() (const Field& F,
								  std::vector<typename Field::Element>& X,
								  const BlasMatrix<typename Field::Element>& A,
								  const std::vector<typename Field::Element>& B) const{
			LQUPMatrix<Field> LQUP(F,A);
			LQUP.right_solve(X,B);
			return X;
		}
		
		std::vector<typename Field::Element>& operator() (const Field& F,
								  const BlasMatrix<typename Field::Element>& A, 
								  std::vector<typename Field::Element>& B) const{
			LQUPMatrix<Field> LQUP(F,A);
			LQUP.right_solve(B);
			return B;
		}
		
	};


	/*
	 * ********************************************************
	 * *** Specialization for TriangularBlasMatrix<Element> ***
	 * ********************************************************
	 */ 


	/*
	 * specialization for Operand of type BlasMatrix<Element>
	 */
		
	template <class Field>
	class BlasMatrixDomainLeftSolve<Field, BlasMatrix<typename Field::Element>,TriangularBlasMatrix<typename Field::Element> > {
	public:
		BlasMatrix<typename Field::Element>& operator() (const Field& F,
								 BlasMatrix<typename Field::Element>& X,
								 const TriangularBlasMatrix<typename Field::Element>& A,
								 const BlasMatrix<typename Field::Element>& B) const{
			
			linbox_check( X.rowdim() == B.rowdim());
			linbox_check( X.coldim() == B.coldim());
			
			typename BlasMatrix<typename Field::Element>::ConstRawIterator  Biter =   B.rawBegin();
			typename BlasMatrix<typename Field::Element>::RawIterator       Xiter =   X.rawBegin();
			
			for (; Biter != B.rawEnd(); ++Biter,++Xiter)
				F.assign(*Xiter,*Biter);		

			return (*this)(F,A,X);
		
	}

		BlasMatrix<typename Field::Element>& operator() (const Field& F,
								 const TriangularBlasMatrix<typename Field::Element>& A, 
								 BlasMatrix<typename Field::Element>& B) const{
			
			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.coldim() == B.rowdim());	
			typename Field::Element _One;
			F.init(_One,1UL);
			
			switch (A.getUpLo()) {
			case BlasTag::up:
				switch(A.getDiag()) {
				case BlasTag::unit:
					FFLAS::ftrsm( F, 
						      FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit,
						      A.rowdim(), B.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
					break;

				case BlasTag::nonunit: {
					//TriangularBlasMatrix<typename Field::Element> Acopy(A);
					FFLAS::ftrsm( F, 
						      FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,
						      A.rowdim(), B.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
					break; }
					
				default:
					throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				}
				break;

			case BlasTag::low:
				switch(A.getDiag()) {
				case BlasTag::unit:
					FFLAS::ftrsm( F, 
						      FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit,
						      A.rowdim(), B.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
					break;

				case BlasTag::nonunit:
					{//TriangularBlasMatrix<typename Field::Element> Acopy(A);
					FFLAS::ftrsm( F, 
						      FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,
						      A.rowdim(), B.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
					break;}

				default:
					throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				}
				break;

			default:
				throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				
			}			
			
			return B;
		}
	};
		


	template <class Field>
	class BlasMatrixDomainRightSolve<Field,BlasMatrix<typename Field::Element>, TriangularBlasMatrix<typename Field::Element> > {
	public:
		BlasMatrix<typename Field::Element>& operator() (const Field& F,
								 BlasMatrix<typename Field::Element>& X,
								 const TriangularBlasMatrix<typename Field::Element>& A,
								 const BlasMatrix<typename Field::Element>& B) const{
			
			linbox_check( X.rowdim() == B.rowdim());
			linbox_check( X.coldim() == B.coldim());

			typename BlasMatrix<typename Field::Element>::ConstRawIterator  Biter =   B.rawBegin();
			typename BlasMatrix<typename Field::Element>::RawIterator       Xiter =   X.rawBegin();
			
			for (; Biter != B.rawEnd(); ++Biter,++Xiter)
				F.assign(*Xiter,*Biter);		
			
			return (*this)(F,A,X);
		}

		BlasMatrix<typename Field::Element>& operator() (const Field& F,
								 const TriangularBlasMatrix<typename Field::Element>& A, 
								 BlasMatrix<typename Field::Element>& B) const{
			
			linbox_check( A.rowdim() == A.coldim());
			linbox_check( B.coldim() == A.rowdim());		
			typename Field::Element _One;
			F.init(_One,1UL);

			switch (A.getUpLo()) {
			case BlasTag::up:
				switch(A.getDiag()) {
				case BlasTag::unit:
					FFLAS::ftrsm( F, 
						      FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit,
						      B.rowdim(), A.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
					break;
				case BlasTag::nonunit:
					FFLAS::ftrsm( F, 
						      FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,
						      B.rowdim(), A.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
					break;
				default:
					throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				}
				break;
			case BlasTag::low:
				switch(A.getDiag()) {
				case BlasTag::unit:
					FFLAS::ftrsm( F, 
						      FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit,
						      B.rowdim(), A.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
					break;
				case BlasTag::nonunit:
					FFLAS::ftrsm( F, 
						      FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,
						      B.rowdim(), A.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
					break;
				default:
					throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				}
				break;
			default:
				throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				
			}
			return B;
		}
	};



	/*
	 * specialization for Operand of type std::vector<Element>
	 */

	template <class Field>
	class BlasMatrixDomainLeftSolve<Field, std::vector<typename Field::Element>, TriangularBlasMatrix<typename Field::Element> > {
	public:
		std::vector<typename Field::Element>& operator() (const Field& F,
								  std::vector<typename Field::Element>& x,
								  const TriangularBlasMatrix<typename Field::Element>& A, 
								  const std::vector<typename Field::Element>& b) const{
			
			linbox_check (x.size() == b.size());
			std::vector<typename Field::Element>::const_iterator biter = b.begin();
			std::vector<typename Field::Element>::iterator       xiter = x.begin();   
			for (;biter!=b.end();++biter,++xiter)
				F.assign(*xiter,*biter);
			
			return (*this)(F,A,x);
		}

		std::vector<typename Field::Element>& operator() (const Field& F,
								  const TriangularBlasMatrix<typename Field::Element>& A,
								  std::vector<typename Field::Element>& b) const{
			
			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.rowdim() == b.size());
					
			switch (A.getUpLo()) {
			case BlasTag::up:
				switch(A.getDiag()) {
				case BlasTag::unit:
					FFLAS::ftrsv( F, 
						      FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit,
						      b.size(),A.getPointer(),A.getStride(),&b[0],1);
					break;
				case BlasTag::nonunit:
					FFLAS::ftrsv( F, 
						      FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,
						      b.size(),A.getPointer(),A.getStride(),&b[0],1);
					break;	
				default:
					throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				}
				break;
			case BlasTag::low:
				switch(A.getDiag()) {
				case BlasTag::unit:
					FFLAS::ftrsv( F, 
						      FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit,
						      b.size(),A.getPointer(),A.getStride(),&b[0],1);
					break;
				case BlasTag::nonunit:
					FFLAS::ftrsv( F, 
						      FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,
						      b.size(),A.getPointer(),A.getStride(),&b[0],1);
					break;
				default:
					throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				}
				break;
			default:
				throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				
			}
			return b;
		}
	};
	
	template <class Field>
	class BlasMatrixDomainRightSolve<Field, std::vector<typename Field::Element>, TriangularBlasMatrix<typename Field::Element> > {
	public:
		std::vector<typename Field::Element>& operator() (const Field& F,
								  std::vector<typename Field::Element>& x,
								  const TriangularBlasMatrix<typename Field::Element>& A,
								  const std::vector<typename Field::Element>& b) const{
			
			linbox_check (x.size() == b.size());
			std::vector<typename Field::Element>::const_iterator biter = b.begin();
			std::vector<typename Field::Element>::iterator       xiter = x.begin();   
			for (;biter!=b.end();++biter,++xiter)
				F.assign(*xiter,*biter);
			
			return (*this)(F,A,x);
		}
		
		std::vector<typename Field::Element>& operator() (const Field& F,
								  const TriangularBlasMatrix<typename Field::Element>& A, 
								  std::vector<typename Field::Element>& b) const{
			
			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.coldim() == b.size());
			
			
			switch (A.getUpLo()) {
			case BlasTag::up:
				switch(A.getDiag()) {
				case BlasTag::unit:
					FFLAS::ftrsv( F, 
						      FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit,
						      b.size(),A.getPointer(),A.getStride(),&b[0],1);
					break;	
				case BlasTag::nonunit:
					FFLAS::ftrsv( F, 
						      FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit,
						      b.size(),A.getPointer(),A.getStride(),&b[0],1);
					break;
				default:
					throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				}
				break;
			case BlasTag::low:
				switch(A.getDiag()) {
				case BlasTag::unit:
					FFLAS::ftrsv( F, 
						      FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit,
						      b.size(),A.getPointer(),A.getStride(),&b[0],1);
					break;	
				case BlasTag::nonunit:
					FFLAS::ftrsv( F, 
						      FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit,
						      b.size(),A.getPointer(),A.getStride(),&b[0],1);
					break;
				default:
					throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				}
				break;
			default:
				throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				
			}
			return b;
		}
	};
	       
	/*
	 *  Method to apply Permutation
	 */
	// Apply a BlasPermutation matrix P to a dense matrix A: 
	// B = A.P 
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyRight(  Operand& B, const Operand& A, const BlasPermutation& P){}

	// B = A.P^t
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyRightTranspose(  Operand& B, const Operand& A, const BlasPermutation& P){}

	// B = P.A 
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyLeft(  Operand& B, const Operand& A, const BlasPermutation& P){}
		
	// B = A.P^t
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyLeftTranspose(  Operand& B, const Operand& A, const BlasPermutation& P){}
		
	// In place apply.
	// A = A.P 
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyinRight( Operand& A, const BlasPermutation& P){}
		
	// A = A.P^t
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyinRightTranspose( Operand& A, const BlasPermutation& P){}       

	// A = P.A 
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyinLeft( Operand& A, const BlasPermutation& P){}
		
	// A = A.P^t
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyinLeftTranspose( Operand& A, const BlasPermutation& P){}

	// Conversion from BlasPermutation to BlackBoxPermutation 
	//template <class Field>
	//inline Permutation& BlasMatrixDomain<Field>::convert ( Permutation& P, const BlasPermutation& BP ){}


} //end of namespace LinBox
	
#endif
