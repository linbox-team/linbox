/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/blas-domain.h
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

#ifndef __BLAS_MATRIX_DOMAIN_H
#define __BLAS_MATRIX_DOMAIN_H

#include <iostream>
#include <vector>
#include <linbox/fflapack/fflapack.h>
#include <linbox/fflas/fflas.h>
#include <linbox/blackbox/permutation.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/util/debug.h>





namespace LinBox {
		 
	/*  Class handling multiplication of a Matrix by an Operand with accumulation an scaling
	 *  Operand can be either a matrix or a vector 
	 *  
	 *  only  function:  operator () are defined :
	 *       D = beta.C + alpha. A*B 
	 *       C = beta.C + alpha. A*B  
	 */
	template< class Field, class Operand1, class Operand2, class Operand3>
	class BlasMatrixDomainMulAdd {
	public:
		Operand1 &operator() (const Field &F, 
				     Operand1 &D,
				     const typename Field::Element &beta, const Operand1 &C,
				     const typename Field::Element &alpha, const Operand2 &A, const Operand3 &B) const;

		Operand1 &operator() (const Field &F,
				     const typename Field::Element &beta, const Operand1 &C,
				     const typename Field::Element &alpha, const Operand2 &A, const Operand3 &B) const;


		// allowing disymetry of Operand2 and Operand3 (only if different type)
		Operand1 &operator() (const Field &F, 
				     Operand1 &D,
				     const typename Field::Element &beta, const Operand1 &C,
				     const typename Field::Element &alpha, const Operand3 &A, const Operand2 &B) const;

		Operand1 &operator() (const Field &F,
				     const typename Field::Element &beta, const Operand1 &C,
				     const typename Field::Element &alpha, const Operand3 &A, const Operand2 &B) const;


		
	};
	
	/*  Class handling in-place multiplication of a Matrix by an Operand 
	 *  Operand can be either a matrix a permutation or a vector 
	 *  
	 *  only  function:  operator () are defined :
	 *       A = A*B
	 *       B = A*B
	 * Note that in-place multiplications are proposed for the specialization 
	 * with a matrix and a permutation.
	 * Using mulin with two matrices is still defined but is non-sense
	 */
	template< class Field, class Operand1, class Operand2, class Operand3>
	class BlasMatrixDomainMul {
	public:
		Operand1 &operator() (const Field &F, 
				      Operand1 &C, const Operand2 &A, const Operand3 &B) const
		{
			typename Field::Element zero, one;
			F.init( zero, 0UL );
			F.init( one, 1UL );
			return BlasMatrixDomainMulAdd<Field,Operand1,Operand2,Operand3>()( F, zero, C, one, A, B );
		}
		
	};
	template< class Field, class Operand1, class Operand2>
	class BlasMatrixDomainMulin {
	public:
		// Defines a dummy mulin over generic matrices using a temporary
		Operand1 &operator() (const Field &F,
				      Operand1 &A, const Operand2 &B) const;
// 		{
// 			typename Field::Element zero, one;
// 			F.init( zero, 0UL );
// 			F.init( one, 1UL );
// 			Operand1* tmp = new Operand1(A);
// 			// Effective copy of A
// 			tmp = A;
// 			BlasMatrixDomainMulAdd<Field,Operand1,Operand1,Operand2>()( F, zero, A, one, *tmp, B );
// 			delete tmp;
// 			return A;
// 		}
		
		Operand2 &operator() (const Field &F, 
				      const Operand2 &A, Operand1 &B ) const;
// 		{
// 			typename Field::Element zero, one;
// 			F.init( zero, 0UL );
// 			F.init( one, 1UL );
// 			Operand2* tmp = new Operand2(B);
// 			// Effective copy of B
// 			*tmp = B;
// 			BlasMatrixDomainMulAdd<Field,Operand2,Operand1,Operand2>()( F, zero, B, one, A, *tmp );
// 			delete tmp;
// 			return B;
// 			}

// 		// allowing disymetry of Operand2 and Operand3 (only if different type)
// 		Operand2 &operator() (const Field &F,
// 				      Operand2 &A, const Operand1 &B) const{
// 			typename Field::Element zero, one;
// 			F.init( zero, 0UL );
// 			F.init( one, 1UL );
// 			Operand2* tmp = new Operand2(A);
// 			// Effective copy of A
// 			tmp = A;
// 			BlasMatrixDomainMulAdd<Field,Operand1,Operand2,Operand1>()( F, zero, A, one, *tmp, B );
// 			delete tmp;
// 			return A;
// 		}
		
// 		Operand1 &operator() (const Field &F, 
// 				      const Operand2 &A, Operand1 &B ) const{
// 			typename Field::Element zero, one;
// 			F.init( zero, 0UL );
// 			F.init( one, 1UL );
// 			Operand1* tmp = new Operand1(B);
// 			// Effective copy of B
// 			*tmp = B;
// 			BlasMatrixDomainMulAdd<Field,Operand1,Operand2,Operand1>()( F, zero, B, one, A, *tmp );
// 			delete tmp;
// 			return B;
// 		}

	};
	
	/*  Class handling inversion of a Matrix 
	 *  
	 *  only  function:  operator () are defined :
	 *       Ainv = A^(-1)
	 *
	 *  Beware, if A is not const this allows an inplace computation
	 *  and so A will be modified
	 */
	template< class Field, class Matrix>
	class BlasMatrixDomainInv {
	public:
		Matrix &operator() (const Field &F, Matrix &Ainv, const Matrix &A) const;
		Matrix &operator() (const Field &F, Matrix &Ainv, Matrix &A) const;
	};
	
	
	
	/*  Class handling rank computation of a Matrix 
	 *  
	 *  only  function:  operator () are defined :
	 *       return the rank of A 
	 * 
	 *  Beware, if A is not const this allows an inplace computation
	 *  and so A will be modified
	 */
	template< class Field, class Matrix>
	class BlasMatrixDomainRank {
	public:
		unsigned int &operator() (const Field &F, const Matrix& A) const;
		unsigned int &operator() (const Field &F, Matrix& A) const;
	};

	/*  Class handling determinant computation of a Matrix 
	 *  
	 *  only  function:  operator () are defined :
	 *       return the determinant of A 
	 * 
	 *  Beware, if A is not const this allows an inplace computation
	 *  and so A will be modified
	 */
	template< class Field, class Matrix>
	class BlasMatrixDomainDet {
	public:
		typename Field::Element operator() (const Field &F, const Matrix& A) const;
		typename Field::Element operator() (const Field &F, Matrix& A) const;
	};


	/*  Class handling resolution of linear system of a Matrix 
	 *  with Operand as right or left and side
	 *  
	 *  only  function:  operator () are defined :
	 *      X = A^(-1).B
	 *      B = A^(-1).B
	 */
	template< class Field, class Operand, class Matrix>
	class BlasMatrixDomainLeftSolve {
	public:
		Operand &operator() (const Field &F, Operand &X, const Matrix &A, const Operand &B) const;
		Operand &operator() (const Field &F, const Matrix &A, Operand &B) const;
	};

	/*  Class handling resolution of linear system of a Matrix 
	 *  with Operand as right or left and side
	 *  
	 *  only  function:  operator () are defined :
	 *      X = B.A^(-1)
	 *      B = B.A^(-1)
	 */
	template< class Field, class Operand, class Matrix>
	class BlasMatrixDomainRightSolve {
	public:
		Operand &operator() (const Field &F, Operand &X, const Matrix &A, const Operand &B) const;
		Operand &operator() (const Field &F, const Matrix &A, Operand &B) const;
	};
	


	/* 
	 *  Interface for all functionnalities provided 
	 *  for BlasMatrix through specialization of all  
	 *  classes defined above.
	 */


	template <class Field>
	class BlasMatrixDomain {

	public:
		typedef typename Field::Element         Element;
		
	private:
    
		Field  _F;
		Element _One;
		Element _Zero;
		Element _MOne;

	public:

		// Constructor of LSP.
		// Initialize the matrix L,S and the permutation P and some constant
		BlasMatrixDomain (const Field& F ) : _F(F) { F.init(_One,1UL); F.init(_Zero,0UL);F.init(_MOne,-1L);}
	    
		// Copy constructor
		BlasMatrixDomain (const BlasMatrixDomain<Field> & BMD) : _F(BMD._F), _One(BMD._One), _Zero(BMD._Zero), _MOne(BMD._MOne) {}


		// Field accessor
		Field& field() {return _F;}

			
		/*
		 * Basics operation available matrix respecting BlasMatrix interface
		 */
 
		// multiplication
		// C = A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& mul(Operand1& C, const Operand2& A, const Operand3& B) const { return BlasMatrixDomainMul<Field,Operand1,Operand2,Operand3>()(_F,C,A,B);}

		

		// multiplication with scaling
		// C = alpha.A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& mul(Operand1& C, const Element& alpha, const Operand2& A, const Operand3& B) const {return muladdin(_Zero,C,alpha,A,B);}

	
		// In place multiplication
		// A = A*B 
		template <class Operand1, class Operand2>
		Operand1& mulin_left(Operand1& A, const Operand2& B ) const { return BlasMatrixDomainMulin<Field,Operand1,Operand2>()(_F,A,B);}
		
		// In place multiplication
		// B = A*B 
		template <class Operand1, class Operand2>
		Operand2& mulin_right(const Operand1& A, Operand2& B ) const { return BlasMatrixDomainMulin<Field,Operand2,Operand1>()(_F,A,B);}

		///////////////////////////////////////////////////////////////////////////////////////////////////
		
// 		// Apply a permutation on the columns of A
// 		// B = A*P
// 		template<class Element >
// 		BlasMatrix<Element>& mul( BlasMatrix<Element>& B, const BlasMatrix<Element>& A, const BlasPermutation& P ){
// 			B = A;
// 			return mul( B, P);
// 		}
		
// 		// Apply a transposed permutation on the columns of A
// 		// B = A*P^t
// 		template<class Element >
// 		BlasMatrix<Element>& mul( BlasMatrix<Element>& B, const BlasMatrix<Element>& A, const TransposedBlasMatrix<BlasPermutation>& P ){
// 			B = A;
// 			return mul( B, P);
// 		}
		
// 		// Apply a permutation on the rows of A
// 		// B = P*A
// 		template<class Element >
// 		BlasMatrix<Element>& mul( BlasMatrix<Element>& B, const BlasPermutation& P, const BlasMatrix<Element>& A ){
// 			B = A;
// 			return mul( P, B );
// 		}
		
// 		// Apply a transposed permutation on the rows of A
// 		// B = A*P^t
// 		template<class Element >
// 		BlasMatrix<Element>& mul( BlasMatrix<Element>& B, const TransposedBlasMatrix<BlasPermutation>& P, const BlasMatrix<Element>& A ){
// 			B = A;
// 			return mul( P, B );
// 		}

// 		// Apply a permutation on the columns of A in place
// 		// A = A*P
// 		template<class Element >
// 		BlasMatrix<Element>& mul( BlasMatrix<Element>& A, const BlasPermutation& P ){
// 			FFLAPACK::applyP( _F, FFLAS::FflasRight, FFLAS::FflasNoTrans, A.rowdim(), 0, A.coldim(), A.getPointer(), A.getStride(), P.getPointer() );
// 			return A;
// 		}
		
// 		// Apply a transposed permutation on the columns of A in place
// 		// A = A*P^t
// 		template<class Element >
// 		BlasMatrix<Element>& mul( BlasMatrix<Element>& A, const TransposedBlasMatrix<BlasPermutation>& P ){
// 			FFLAPACK::applyP( _F, FFLAS::FflasRight, FFLAS::FflasTrans, A.rowdim(), 0, A.coldim(), A.getPointer(), A.getStride(), P.getMatrix().getPointer() );
// 			return A;
// 		}
// 		// Apply a permutation on the rows of A in place
// 		// B = P*A
// 		template<class Element >
// 		BlasMatrix<Element>& mul( const BlasPermutation& P, BlasMatrix<Element>& A ){
// 			FFLAPACK::applyP( _F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, A.coldim(), 0, A.rowdim(), A.getPointer(), A.getStride(), P.getPointer() );
//  			return A;
// 		}
		
// 		// Apply a transposed permutation on the rows of A in place
// 		// A = P^t*A
// 		template<class Element >
// 		BlasMatrix<Element>& mul( const TransposedBlasMatrix<BlasPermutation>& P, BlasMatrix<Element>& A ){
// 			FFLAPACK::applyP( _F, FFLAS::FflasLeft, FFLAS::FflasTrans, A.coldim(), 0, A.rowdim(), A.getPointer(), A.getStride(), P.getMatrix().getPointer() );
// 			return A;
// 		}

		///////////////////////////////////////////////////////////////////////////////////////////////////
		// axpy
		// D = A*B + C
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axpy(Operand1& D, const Operand2& A, const Operand3& B, const Operand1& C) const {return muladd(D,_One,C,_One,A,B);}

		// axpyin
		// C += A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axpyin(Operand1& C, const Operand2& A, const Operand3& B) const {return muladdin(_One,C,_One,A,B);}
 
		// axmy
		// D= A*B - C
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axmy(Operand1& D, const Operand2& A, const Operand3& B, const Operand1& C) const {return muladd(D,_MOne,C,_One,A,B);}

		// axmyin
		// C = A*B - C
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axmyin(Operand1& C, const Operand2& A, const Operand3& B) const {return muladdin(_MOne,C,_One,A,B);}
		
		//  general matrix-matrix multiplication and addition with scaling
		// D= beta.C + alpha.A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& muladd(Operand1& D, const Element& beta, const Operand1& C,
				const Element& alpha, const Operand2& A, const Operand3& B) const {
			return BlasMatrixDomainMulAdd<Field,Operand1,Operand2,Operand3>()(_F,D,beta,C,alpha,A,B);
		}
		
		// C= beta.C + alpha.A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& muladdin(const Element& beta, Operand1& C,
				   const Element& alpha, const Operand2& A, const Operand3& B) const {
			return BlasMatrixDomainMulAdd<Field,Operand1,Operand2,Operand3>()(_F,beta,C,alpha,A,B);
		}


		/*
		 * Solutions available for matrix respecting BlasMatrix interface
		 */	

		// Inversion
		template <class Matrix>
		Matrix& inv( Matrix &Ainv, const Matrix &A) const { 
			return BlasMatrixDomainInv<Field,Matrix>()(_F,Ainv,A);
		}
		
		// Inversion (the matrix A is modified)
		template <class Matrix>
		Matrix& invin( Matrix &Ainv, Matrix &A) const { 
			return BlasMatrixDomainInv<Field,Matrix>()(_F,Ainv,A);
		}

		// Rank
		template <class Matrix>
		unsigned int rank(const Matrix &A) const {
			return BlasMatrixDomainRank<Field,Matrix>()(_F,A);
		}

		// in-place Rank (the matrix is modified)
		template <class Matrix>
		unsigned int rankin(Matrix &A) const {
			return BlasMatrixDomainRank<Field, Matrix>()(_F,A);
		}

		// determinant
		template <class Matrix>
		Element det(const Matrix &A) const {
			return BlasMatrixDomainDet<Field, Matrix>()(_F,A);
		}

		//in-place Determinant (the matrix is modified)
		template <class Matrix>
		Element detin(Matrix &A) const {
			return BlasMatrixDomainDet<Field, Matrix>()(_F,A);
		}
		
		/*
		 * Solvers for Matrix (respecting BlasMatrix interface) 
		 * with right or left Operand hand side
		 */ 

		// non-singular linear solve with matrix right hand side 
		// AX=B
		template <class Operand, class Matrix>
		Operand& left_solve (Operand& X, const Matrix& A, const Operand& B) const {
			return BlasMatrixDomainLeftSolve<Field,Operand,Matrix>()(_F,X,A,B);
		}
		
		// non-singular linear solve with matrix right hand side, the result is stored in-place in B
		// AX=B , (B<-X)
		template <class Operand,class Matrix>
		Operand& left_solve (const Matrix& A, Operand& B) const {
			return BlasMatrixDomainLeftSolve<Field,Operand,Matrix>()(_F,A,B);
		}
		
		// non-singular linear solve with matrix right hand side 
		// XA=B
		template <class Operand, class Matrix>
		Operand& right_solve (Operand& X, const Matrix& A, const Operand& B) const {
			return BlasMatrixDomainRightSolve<Field,Operand,Matrix>()(_F,X,A,B);
		}
		
		// non-singular linear solve with matrix right hand side, the result is stored in-place in B
		// XA=B , (B<-X)
		template <class Operand, class Matrix>
		Operand& right_solve (const Matrix& A, Operand& B) const {
			return BlasMatrixDomainRightSolve<Field,Operand,Matrix>()(_F,A,B);
		}
		
	}; /* end of class BlasMatrixDomain */


} /* end of namespace LinBox */

#include <linbox/algorithms/blas-domain.inl>

#endif /* __BLAS_DOMAIN_H*/

