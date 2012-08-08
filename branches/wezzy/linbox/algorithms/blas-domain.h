/* linbox/algorithms/blas-domain.h
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
 *               Clément Pernet clement.pernet@imag.fr
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file algorithms/blas-domain.h
 * @ingroup algorithms
 * @brief NO DOC
 * @warning A <code>BlasMatrixDomain<Field></code> should be templated by a
 * \link LinBox::Modular Modular\endlink field. In particular, this domain
 * is not suitable for integers.
 * @warning A \e Field does mean here a \e Field and not a general \f$\mathbf{Z}/m\mathbf{Z}\f$ \e ring. You'll be warned...
 */

#ifndef __LINBOX_blas_matrix_domain_H
#define __LINBOX_blas_matrix_domain_H

#include <iostream>
#include <vector>

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#include <fflas-ffpack/ffpack/ffpack.h>
#include <fflas-ffpack/fflas/fflas.h>
#include "linbox/matrix/blas-matrix.h"
#include "linbox/matrix/permutation-matrix.h"


namespace LinBox
{

	const int BlasBound = 1 << 26;

	/** @internal
	 * Class handling multiplication of a Matrix by an Operand with accumulation and scaling.
	 *  Operand can be either a matrix or a vector.
	 *
	 *  The only function:  operator () is defined :
	 *     -  D = beta.C + alpha. A*B
	 *     -  C = beta.C + alpha. A*B
	 */
	template< class Field, class Operand1, class Operand2, class Operand3>
	class BlasMatrixDomainMulAdd {
	public:
		Operand1 &operator() (const Field &F,
				      Operand1 &D,
				      const typename Field::Element &beta, const Operand1 &C,
				      const typename Field::Element &alpha, const Operand2 &A, const Operand3 &B) const;

		Operand1 &operator() (const Field &F,
				      const typename Field::Element &beta, Operand1 &C,
				      const typename Field::Element &alpha, const Operand2 &A, const Operand3 &B) const;


		// allowing disymetry of Operand2 and Operand3 (only if different type)
		Operand1 &operator() (const Field &F,
				      Operand1 &D,
				      const typename Field::Element &beta, const Operand1 &C,
				      const typename Field::Element &alpha, const Operand3 &A, const Operand2 &B) const;

		Operand1 &operator() (const Field &F,
				      const typename Field::Element &beta, Operand1 &C,
				      const typename Field::Element &alpha, const Operand3 &A, const Operand2 &B) const;



	};

	/*!@internal
	 * Adding two matrices
	 */
	template< class Field, class Operand1, class Operand2, class Operand3>
	class BlasMatrixDomainAdd {
	public:
		Operand1 &operator() (const Field &F,
				      Operand1 &C, const Operand2 &A, const Operand3 &B) const;

	};

	/*!@internal
	 * Substracting two matrices
	 */
	template< class Field, class Operand1, class Operand2, class Operand3>
	class BlasMatrixDomainSub {
	public:
		Operand1 &operator() (const Field &F,
				      Operand1 &C, const Operand2 &A, const Operand3 &B) const;

	};
	//! C -= A
	template< class Field, class Operand1, class Operand2>
	class BlasMatrixDomainSubin {
	public:
		Operand1 &operator() (const Field &F,
				      Operand1 &C, const Operand2 &A) const;

	};
	//! C += A
	template< class Field, class Operand1, class Operand2>
	class BlasMatrixDomainAddin {
	public:
		Operand1 &operator() (const Field &F,
				      Operand1 &C, const Operand2 &A) const;

	};


	/*!@internal
	 * Copying a matrix
	 */
	template< class Field, class Operand1, class Operand2>
	class BlasMatrixDomainCopy {
	public:
		Operand1 &operator() (const Field &F,
				      Operand1 &B, const Operand2 &A) const;

	};

	/*!@internal
	 * multiplying two matrices.
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

	/*! @internal
	 * Class handling in-place multiplication of a Matrix by an Operand.
	 *  Operand can be either a matrix a permutation or a vector
	 *
	 *  only  function:  operator () are defined :
	 *     -  A = A*B
	 *     -  B = A*B
	 *     .
	 * @note In-place multiplications are proposed for the specialization
	 * with a matrix and a permutation.
	 * @warning Using mulin with two matrices is still defined but is
	 * non-sense
	 */
	// Operand 2 is always the type of the matrix which is not modified
	// ( for example: BlasPermutation TriangularBlasMatrix )
	template< class Field, class Operand1, class Operand2>
	class BlasMatrixDomainMulin {
	public:
		// Defines a dummy mulin over generic matrices using a temporary
		Operand1 &operator() (const Field &F,
				      Operand1 &A, const Operand2 &B) const
		{
			typename Field::Element zero, one;
			F.init( zero, 0UL );
			F.init( one, 1UL );
			Operand1* tmp = new Operand1(A);
			// Effective copy of A
			*tmp = A;
			BlasMatrixDomainMulAdd<Field,Operand1,Operand1,Operand2>()( F, zero, A, one, *tmp, B );
			delete tmp;
			return A;
		}

		Operand1 &operator() (const Field &F,
				      const Operand2 &A, Operand1 &B ) const
		{
			typename Field::Element zero, one;
			F.init( zero, 0UL );
			F.init( one, 1UL );
			Operand1* tmp = new Operand1(B);
			// Effective copy of B
			*tmp = B;
			BlasMatrixDomainMulAdd<Field,Operand1,Operand1,Operand2>()( F, zero, B, one, A, *tmp );
			delete tmp;
			return B;
		}
	};

	/*! @internal
	 * Class handling inversion of a Matrix.
	 *
	 *  only  function:  operator () are defined :
	 *    -   Ainv = A^(-1)
	 *
	 *  Returns nullity of matrix (0 iff inversion was ok)
	 *
	 * @warning Beware, if A is not const this allows an inplace computation
	 *  and so A will be modified
	 *
	 */
	template< class Field, class Matrix>
	class BlasMatrixDomainInv {
	public:
		int &operator() (const Field &F, Matrix &Ainv, const Matrix &A) const;
		int &operator() (const Field &F, Matrix &Ainv, Matrix &A) const;
	};

	/*! @internal
	 * Class handling rank computation of a Matrix.
	 *
	 *  only  function:  operator () are defined :
	 *    -  return the rank of A
	 *
	 *  @warning Beware, if A is not const this allows an inplace computation
	 *  and so A will be modified
	 */
	template< class Field, class Matrix>
	class BlasMatrixDomainRank {
	public:
		unsigned int &operator() (const Field &F, const Matrix& A) const;
		unsigned int &operator() (const Field &F, Matrix& A) const;
	};

	/*! @internal
	 * Class handling determinant computation of a Matrix.
	 *
	 *  only  function:  operator () are defined :
	 *     -  return the determinant of A
	 *
	 * @warning Beware, if A is not const this allows an inplace computation
	 *  and so A will be modified
	 */
	template< class Field, class Matrix>
	class BlasMatrixDomainDet {
	public:
		typename Field::Element operator() (const Field &F, const Matrix& A) const;
		typename Field::Element operator() (const Field &F, Matrix& A) const;
	};

	/*! @internal
	 * Class handling resolution of linear system of a Matrix.
	 *  with Operand as right or left and side
	 *
	 *  only  function:  operator () are defined :
	 *    -  X = A^(-1).B
	 *    -  B = A^(-1).B
	 */
	template< class Field, class Operand, class Matrix>
	class BlasMatrixDomainLeftSolve {
	public:
		Operand &operator() (const Field &F, Operand &X, const Matrix &A, const Operand &B) const;
		Operand &operator() (const Field &F, const Matrix &A, Operand &B) const;
	};

	/*! @internal
	 * Class handling resolution of linear system of a Matrix.
	 *  with Operand as right or left and side
	 *
	 *  only  function:  operator () are defined :
	 *   -   X = B.A^(-1)
	 *   -   B = B.A^(-1)
	 */
	template< class Field, class Operand, class Matrix>
	class BlasMatrixDomainRightSolve {
	public:
		Operand &operator() (const Field &F, Operand &X, const Matrix &A, const Operand &B) const;
		Operand &operator() (const Field &F, const Matrix &A, Operand &B) const;
	};

	/*! @internal
	 * Class handling the minimal polynomial  of a Matrix.
	 */
	template< class Field, class Polynomial, class Matrix>
	class BlasMatrixDomainMinpoly {
	public:
		Polynomial&  operator() (const Field &F, Polynomial& P, const Matrix& A) const;
	};

	/*! @internal
	 * Class handling the characteristic polynomial  of a Matrix.
	 */
	template< class Field, class ContPol, class Matrix>
	class BlasMatrixDomainCharpoly {
	public:
		//		typedef Container<Polynomial> ContPol;

		ContPol&  operator() (const Field &F, ContPol& P, const Matrix& A) const;
	};

	/**
	 *  Interface for all functionnalities provided
	 *  for BlasMatrix.
	 *  @internal
	 *  Done through specialization of all
	 *  classes defined above.
	 */
	template <class Field_>
	class BlasMatrixDomain {

	public:
		typedef Field_ Field;
		typedef typename Field::Element         Element;
		typedef BlasSubmatrix<Field> Matrix;

	protected:

		const Field  * _field;
		Element _One;
		Element _Zero;
		Element _MOne;

	public:

		//! Constructor of BlasDomain.

		BlasMatrixDomain (const Field& F ) :
			_field(&F)
		{
			F.init(_One,1UL);
			F.init(_Zero,0UL);
			F.init(_MOne,-1);
#ifndef NDEBUG
			if (!Givaro::probab_prime(F.characteristic())) {
				std::cout << " *** WARNING *** "                                           << std::endl;
				std::cout << " You are using a BLAS Domain where your field is not prime " << 
				F.characteristic() << std::endl;
			}
#endif
		}

		//! Copy constructor
		BlasMatrixDomain (const BlasMatrixDomain<Field> & BMD) :
			_field(BMD._field), _One(BMD._One), _Zero(BMD._Zero), _MOne(BMD._MOne)
		{
#ifndef NDEBUG
			if (!Givaro::probab_prime(field().characteristic())) {
				std::cout << " *** WARNING *** "                                           << std::endl;
				std::cout << " You are using a BLAS Domain where your field is not prime " << std::endl;
			}
#endif

		}


		//! Field accessor
		const Field& field() const { return *_field; }

		/*
		 * Basics operation available matrix respecting BlasMatrix interface
		 */

		//! multiplication.
		//! C = A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& mul(Operand1& C, const Operand2& A, const Operand3& B) const
		{
			return BlasMatrixDomainMul<Field,Operand1,Operand2,Operand3>()(field(),C,A,B);
		}

		//! addition.
		//! C = A+B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& add(Operand1& C, const Operand2& A, const Operand3& B) const
		{
			return BlasMatrixDomainAdd<Field,Operand1,Operand2,Operand3>()(field(),C,A,B);
		}

		//! copy.
		//! B = A
		template <class Operand1, class Operand2>
		Operand1& copy(Operand1& B, const Operand2& A) const
		{
			return BlasMatrixDomainCopy<Field,Operand1,Operand2>()(field(),B,A);
		}

		//! substraction
		//! C = A-B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& sub(Operand1& C, const Operand2& A, const Operand3& B) const
		{
			return BlasMatrixDomainSub<Field,Operand1,Operand2,Operand3>()(field(),C,A,B);
		}

		//! substraction (in place)
		//! C -= B
		template <class Operand1, class Operand3>
		Operand1& subin(Operand1& C, const Operand3& B) const
		{
			return BlasMatrixDomainSubin<Field,Operand1,Operand3>()(field(),C,B);
		}

		//! addition (in place)
		//! C += B
		template <class Operand1, class Operand3>
		Operand1& addin(Operand1& C, const Operand3& B) const
		{
			return BlasMatrixDomainAddin<Field,Operand1,Operand3>()(field(),C,B);
		}


		//! multiplication with scaling.
		//! C = alpha.A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& mul(Operand1& C, const Element& alpha, const Operand2& A, const Operand3& B) const
		{
			return muladdin(_Zero,C,alpha,A,B);
		}


		//! In place multiplication.
		//! A = A*B
		template <class Operand1, class Operand2>
		Operand1& mulin_left(Operand1& A, const Operand2& B ) const
		{
			return BlasMatrixDomainMulin<Field,Operand1,Operand2>()(field(),A,B);
		}

		//! In place multiplication.
		//! B = A*B
		template <class Operand1, class Operand2>
		Operand2& mulin_right(const Operand1& A, Operand2& B ) const
		{
			return BlasMatrixDomainMulin<Field,Operand2,Operand1>()(field(),A,B);
		}

		//! axpy.
		//! D = A*B + C
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axpy(Operand1& D, const Operand2& A, const Operand3& B, const Operand1& C) const
		{
			return muladd(D,_One,C,_One,A,B);
		}

		//! axpyin.
		//! C += A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axpyin(Operand1& C, const Operand2& A, const Operand3& B) const
		{
			return muladdin(_One,C,_One,A,B);
		}

		//! maxpy.
		//! D = C - A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& maxpy(Operand1& D, const Operand2& A, const Operand3& B, const Operand1& C)const
		{
			return muladd(D,_One,C,_MOne,A,B);
		}

		//! maxpyin.
		//! C -= A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& maxpyin(Operand1& C, const Operand2& A, const Operand3& B) const
		{
			return muladdin(_One,C,_MOne,A,B);
		}

		//! axmy.
		//! D= A*B - C
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axmy(Operand1& D, const Operand2& A, const Operand3& B, const Operand1& C) const
		{
			return muladd(D,_MOne,C,_One,A,B);
		}

		//! axmyin.
		//! C = A*B - C
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axmyin(Operand1& C, const Operand2& A, const Operand3& B) const
		{
			return muladdin(_MOne,C,_One,A,B);
		}

		//!  general matrix-matrix multiplication and addition with scaling.
		//! D= beta.C + alpha.A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& muladd(Operand1& D, const Element& beta, const Operand1& C,
				 const Element& alpha, const Operand2& A, const Operand3& B) const
		{
			return BlasMatrixDomainMulAdd<Field,Operand1,Operand2,Operand3>()(field(),D,beta,C,alpha,A,B);
		}

		//! muladdin.
		//! C= beta.C + alpha.A*B.
		template <class Operand1, class Operand2, class Operand3>
		Operand1& muladdin(const Element& beta, Operand1& C,
				   const Element& alpha, const Operand2& A, const Operand3& B) const
		{
			return BlasMatrixDomainMulAdd<Field,Operand1,Operand2,Operand3>()(field(),beta,C,alpha,A,B);
		}


		/*!
		 * @name Solutions available for matrix respecting BlasMatrix interface
		 */
		//@{

		//! Inversion
		template <class Matrix>
		Matrix& inv( Matrix &Ainv, const Matrix &A) const
		{
			BlasMatrixDomainInv<Field,Matrix>()(field(),Ainv,A);
			return Ainv;
		}

		//! Inversion (in place)
		template <class Matrix>
		Matrix& invin( Matrix &Ainv, Matrix &A) const
		{
			BlasMatrixDomainInv<Field,Matrix>()(field(),Ainv,A);
			return Ainv;
		}

		//! Inversion (the matrix A is modified)
		template <class Matrix>
		Matrix& invin(Matrix &A) const
		{
			Matrix tmp(A.rowdim(), A.coldim());
			tmp = A;
			BlasMatrixDomainInv<Field,Matrix>()(field(),A,tmp);
			return A;
		}


		/*! Division.
		 * C = A B^{-1}  ==>  C . B = A
		 */
		template <class Matrix>
		Matrix& div( Matrix &C, const Matrix &A, const Matrix &B) const
		{
			return this->right_solve(C,B,A);
		}


		//! Inversion w singular check
		template <class Matrix>
		Matrix& inv( Matrix &Ainv, const Matrix &A, int& nullity) const
		{
			nullity = BlasMatrixDomainInv<Field,Matrix>()(field(),Ainv,A);
			return Ainv;
		}

		//! Inversion (the matrix A is modified) w singular check
		template <class Matrix>
		Matrix& invin( Matrix &Ainv, Matrix &A, int& nullity) const
		{
			nullity = BlasMatrixDomainInv<Field,Matrix>()(field(),Ainv,A);
			return Ainv;
		}

		//! Rank
		template <class Matrix>
		unsigned int rank(const Matrix &A) const
		{
			return BlasMatrixDomainRank<Field,Matrix>()(field(),A);
		}

		//! in-place Rank (the matrix is modified)
		template <class Matrix>
		unsigned int rankin(Matrix &A) const
		{
			return BlasMatrixDomainRank<Field, Matrix>()(field(),A);
		}

		//! determinant
		template <class Matrix>
		Element det(const Matrix &A) const
		{
			return BlasMatrixDomainDet<Field, Matrix>()(field(),A);
		}

		//! in-place Determinant (the matrix is modified)
		template <class Matrix>
		Element detin(Matrix &A) const
		{
			return BlasMatrixDomainDet<Field, Matrix>()(field(),A);
		}
		//@}

		/*!
		 * @name Solvers for Matrix (respecting BlasMatrix interface)
		 * with Operand as right or left hand side
		 */
		//@{
		//! linear solve with matrix right hand side.
		//! AX=B
		template <class Operand, class Matrix>
		Operand& left_solve (Operand& X, const Matrix& A, const Operand& B) const
		{
			return BlasMatrixDomainLeftSolve<Field,Operand,Matrix>()(field(),X,A,B);
		}

		//! linear solve with matrix right hand side, the result is stored in-place in B.
		//! @pre A must be square
		//! AX=B , (B<-X)
		template <class Operand,class Matrix>
		Operand& left_solve (const Matrix& A, Operand& B) const
		{
			return BlasMatrixDomainLeftSolve<Field,Operand,Matrix>()(field(),A,B);
		}

		//! linear solve with matrix right hand side.
		//! XA=B
		template <class Operand, class Matrix>
		Operand& right_solve (Operand& X, const Matrix& A, const Operand& B) const
		{
			return BlasMatrixDomainRightSolve<Field,Operand,Matrix>()(field(),X,A,B);
		}

		//! linear solve with matrix right hand side, the result is stored in-place in B.
		//! @pre A must be square
		//! XA=B , (B<-X)
		template <class Operand, class Matrix>
		Operand& right_solve (const Matrix& A, Operand& B) const
		{
			return BlasMatrixDomainRightSolve<Field,Operand,Matrix>()(field(),A,B);
		}

		//! minimal polynomial computation.
		template <class Polynomial, class Matrix>
		Polynomial& minpoly( Polynomial& P, const Matrix& A ) const
		{
			return BlasMatrixDomainMinpoly<Field, Polynomial, Matrix>()(field(),P,A);
		}

		//! characteristic polynomial computation.
		template <class Polynomial,  class Matrix >
		Polynomial& charpoly( Polynomial& P, const Matrix& A ) const
		{

			commentator().start ("Modular Dense Charpoly ", "MDCharpoly");
			std::list<Polynomial> P_list;
			P_list.clear();
			BlasMatrixDomainCharpoly<Field, std::list<Polynomial>, Matrix >()(field(),P_list,A);


			Polynomial tmp(A.rowdim()+1);
			typename std::list<Polynomial>::iterator it = P_list.begin();
			P = *(it++);
			while( it!=P_list.end() ){
				// Waiting for an implementation of a domain of polynomials
				mulpoly( tmp, P, *it);
				P = tmp;
				//	delete &(*it);
				++it;
			}
			commentator().stop ("done", NULL, "MDCharpoly");

			return P;
		}

		//! characteristic polynomial computation.
		template <class Polynomial, class Matrix >
		std::list<Polynomial>& charpoly( std::list<Polynomial>& P, const Matrix& A ) const
		{
			return BlasMatrixDomainCharpoly<Field, std::list<Polynomial>, Matrix >()(field(),P,A);
		}


		//private:
		//! @todo Temporary: waiting for an implementation of a domain of polynomial
		template<class Polynomial>
		Polynomial &
		mulpoly(Polynomial &res, const Polynomial & P1, const Polynomial & P2)const
		{
			size_t i,j;
			res.resize(P1.size()+P2.size()-1);
			for (i=0;i<res.size();i++)
				field().assign(res[i],_Zero);
			for ( i=0;i<P1.size();i++)
				for ( j=0;j<P2.size();j++)
					field().axpyin(res[i+j],P1[i],P2[j]);
			return res;

		}
		//@}

		template<class Matrix1, class Matrix2>
		bool areEqual(const Matrix1 & A, const Matrix2 & B)
		{
			if ( (A.rowdim() != B.rowdim()) || (A.coldim() != B.coldim()) )
				return false ;
			for (size_t i = 0 ; i < A.rowdim() ; ++i)
				for (size_t j = 0 ; j < A.coldim() ; ++j)
					if (!field().areEqual(A.getEntry(i,j),B.getEntry(i,j))) //!@bug use refs
						return false ;
			return true ;
		}

		template<class Matrix>
		void setIdentity(Matrix & I)
		{
			for (size_t i = 0 ; i< I.rowdim() ; ++i)
				for (size_t j = 0 ; j < I.coldim() ; ++j) {
					if (i == j)
						I.setEntry(i,j,_One);
					else
						I.setEntry(i,j,_Zero);
				}

		}

		template<class Matrix>
		void setZero(Matrix & I)
		{
			// use Iterator
			for (size_t i = 0 ; i< I.rowdim() ; ++i)
				for (size_t j = 0 ; j < I.coldim() ; ++j) {
						I.setEntry(i,j,_Zero);
				}
		}


		template<class Matrix1>
		bool isZero(const Matrix1 & A)
		{
			for (size_t i = 0 ; i < A.rowdim() ; ++i)
				for (size_t j = 0 ; j < A.coldim() ; ++j)
					if (!field().isZero(A.getEntry(i,j))) //!@bug use refs
						return false ;
			return true ;
		}

		template<class Matrix1>
		bool isIdentity(const Matrix1 & A)
		{
			if (A.rowdim() != A.coldim())
				return false ;
			for (size_t i = 0 ; i < A.rowdim() ; ++i)
				if (!field().isOne(A.getEntry(i,i)))
					return false;

			for (size_t i = 0 ; i < A.rowdim() ; ++i)
				for (size_t j = 0 ; j < i ; ++j)
					if (!field().isZero(A.getEntry(i,j))) //!@bug use refs
						return false ;
			for (size_t i = 0 ; i < A.rowdim() ; ++i)
				for (size_t j = i+1 ; j < A.coldim() ; ++j)
					if (!field().isZero(A.getEntry(i,j))) //!@bug use refs
						return false ;
			return true ;
		}

		template<class Matrix1>
		bool isIdentityGeneralized(const Matrix1 & A)
		{
			size_t mn = std::min(A.rowdim(),A.coldim());
			for (size_t i = 0 ; i < mn ; ++i)
				if (!field().isOne(A.getEntry(i,i)))
					return false;

			for (size_t i = 0 ; i < A.rowdim() ; ++i)
				for (size_t j = 0 ; j < std::min(i,mn) ; ++j)
					if (!field().isZero(A.getEntry(i,j))) //!@bug use refs
						return false ;
			for (size_t i = 0 ; i < A.rowdim() ; ++i)
				for (size_t j = i+1 ; j < A.coldim() ; ++j)
					if (!field().isZero(A.getEntry(i,j))) //!@bug use refs
						return false ;
			return true ;
		}

	public:

		/** Print matrix.
		 * @param  os  Output stream to which matrix is written.
		 * @param  A   Matrix.
		 * @returns reference to os.
		 */
		template <class Matrix>
		inline std::ostream &write (std::ostream &os, const Matrix &A) const
		{
			return A.write (os);
		}

		template <class Matrix>
		inline std::ostream &write (std::ostream &os, const Matrix &A, bool maple_format) const
		{
			return A.write (os, field(), maple_format);
		}

		/** Read matrix
		 * @param  is  Input stream from which matrix is read.
		 * @param  A   Matrix.
		 * @returns reference to is.
		 */
		template <class Matrix>
		inline std::istream &read (std::istream &is, Matrix &A) const
		{
			return A.read (is, field());
		}

	}; /* end of class BlasMatrixDomain */


} /* end of namespace LinBox */

#include "linbox/algorithms/blas-domain.inl"

#endif /* __LINBOX_blas_matrix_domain_H */


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

