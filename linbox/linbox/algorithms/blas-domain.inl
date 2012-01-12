/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/blas-domain.inl
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



#ifndef __LINBOX_blas_matrix_domain_INL
#define __LINBOX_blas_matrix_domain_INL

#include "linbox/matrix/blas-matrix.h"
#include "linbox/matrix/factorized-matrix.h"

namespace LinBox
{

	/*
	 * **********************************************
	 * *** Specialization for BlasMatrix<Field> ***
	 * **********************************************
	 */


#if 0
	// Inversion
	// dpritcha: now returns nullity. (2004-07-19)
	// previously returned Ainv but this is passed back anyway.
	template <class Field>
	class BlasMatrixDomainInv<Field,BlasMatrix<Field> > {
	public:
		int operator() (const Field                   &F,
				BlasMatrix<Field>        &Ainv,
				const BlasMatrix<Field>     &A) const
		{

			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.rowdim() == Ainv.rowdim());
			linbox_check( A.coldim() == Ainv.coldim());
			BlasMatrix<Field> tmp(A);
			return (*this)(F,Ainv,tmp);
		}

		int operator() (const Field                &F,
				BlasMatrix<Field>     &Ainv,
				BlasMatrix<Field>        &A) const
		{

			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.rowdim() == Ainv.rowdim());
			linbox_check( A.coldim() == Ainv.coldim());
			int nullity;
			FFPACK::Invert(F,A.rowdim(),A.getPointer(),A.getStride(),
				       Ainv.getPointer(),Ainv.getStride(),nullity);
			return nullity;
		}

	};
#endif


#if 0
	// Rank
	template <class Field>
	class 	BlasMatrixDomainRank<Field,BlasMatrix<Field> > {
	public:
		inline unsigned int operator() (const Field                &F,
						const BlasMatrix<Field>  &A) const
		{

			BlasMatrix<Field> tmp(A);
			return (*this)(F,tmp);
		}

		inline unsigned int operator() (const Field                &F,
						BlasMatrix<Field>        &A) const
		{

			return FFPACK::Rank(F, A.rowdim(), A.coldim(),A.getPointer(), A.getStride());
		}
	};

	// determinant
	template <class Field>
	class BlasMatrixDomainDet<Field,BlasMatrix<Field> > {
	public:
		inline typename Field::Element operator()(const Field                 &F,
							  const BlasMatrix<Field>   &A) const
		{

			BlasMatrix<Field> tmp(A);
			return  (*this)(F,tmp);
		}

		inline typename Field::Element operator() (const Field                &F,
							   BlasMatrix<Field>        &A) const
		{

			return FFPACK::Det(F, A.rowdim(), A.coldim(),A.getPointer(), A.getStride());
		}
	};
#endif


	/*
	 * **********************************************
	 * *** Specialization for BlasMatrix<Field> ***
	 * **********************************************
	 */


	// Inversion
	// dpritcha: now returns nullity. (2004-07-19)
	// previously returned Ainv but this is passed back anyway.
	template <class Field>
	class BlasMatrixDomainInv<Field,BlasMatrix<Field> > {
	public:
		int operator() (const Field                                   &F,
				BlasMatrix<Field>        &Ainv,
				const BlasMatrix<Field>     &A) const
		{

			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.rowdim() == Ainv.rowdim());
			linbox_check( A.coldim() == Ainv.coldim());
			BlasMatrix<Field> tmp(A);
			return (*this)(F,Ainv,tmp);
		}

		int operator() (const Field                                  &F,
				BlasMatrix<Field>       &Ainv,
				BlasMatrix<Field>          &A) const
		{

			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.rowdim() == Ainv.rowdim());
			linbox_check( A.coldim() == Ainv.coldim());
			int nullity;
			FFPACK::Invert (F, A.rowdim(), A.getPointer(), A.getStride(),
					Ainv.getPointer(),Ainv.getStride(),nullity);
			return nullity;
		}

	};

	// Rank
	template <class Field>
	class 	BlasMatrixDomainRank<Field,BlasMatrix<Field> > {
	public:
		inline unsigned int operator() (const Field                                &F,
						const BlasMatrix<Field>  &A) const
		{

			BlasMatrix<Field> tmp(A);
			return (*this)(F,tmp);
		}

		inline unsigned int
		operator() (const Field                           &F,
			    BlasMatrix<Field>   &A) const
		{

			return (unsigned int) FFPACK::Rank(F, A.rowdim(), A.coldim(),A.getPointer(), A.getStride());
		}
	};

	// determinant
	template <class Field>
	class BlasMatrixDomainDet<Field,BlasMatrix<Field> > {
	public:
		inline typename Field::Element operator()(const Field                                &F,
							  const BlasMatrix<Field>  &A) const
		{

			BlasMatrix<Field> tmp(A);
			return  (*this)(F,tmp);
		}

		inline typename Field::Element operator() (const Field                             &F,
							   BlasMatrix<Field>     &A) const
		{

			return FFPACK::Det(F, A.rowdim(), A.coldim(),A.getPointer(), A.getStride());
		}
	};


	/*
	 * specialization for Operand1, Operand2 and Operand3  of type BlasMatrix<Field>
	 */

	template<class Field>
	class 	BlasMatrixDomainAdd<Field,BlasMatrix<Field>,BlasMatrix<Field>, BlasMatrix<Field> > {
	public:
		BlasMatrix<Field>& operator()(const Field& F,
								BlasMatrix<Field>& C,
								const BlasMatrix<Field>& A,
								const BlasMatrix<Field>& B) const
		{
			linbox_check( A.rowdim() == B.rowdim());
			linbox_check( C.rowdim() == A.rowdim());
			linbox_check( A.coldim() == B.coldim());
			linbox_check( C.coldim() == A.coldim());
			FFLAS::fadd (F, C.rowdim(), C.coldim(),
				     A.getPointer(), A.getStride(),
				     B.getPointer(), B.getStride(),
				     C.getPointer(), C.getStride());
			return C;
		}
	};

	template<class Field>
	class 	BlasMatrixDomainCopy<Field,BlasMatrix<Field>, BlasMatrix<Field> > {
	public:
		BlasMatrix<Field>& operator()(const Field& F,
								BlasMatrix<Field>& B,
								const BlasMatrix<Field>& A) const
		{
			linbox_check( A.rowdim() == B.rowdim());
			linbox_check( A.coldim() == B.coldim());
			for (size_t i=0; i<A.rowdim(); i++)
				FFLAS::fcopy (F, A.coldim(),
					      B.getPointer() + i*B.getStride(), 1,
					      A.getPointer() + i*A.getStride(), 1);
			return B;
		}
	};

	template<class Field>
	class 	BlasMatrixDomainSub<Field,BlasMatrix<Field>,BlasMatrix<Field>, BlasMatrix<Field> > {
	public:
		BlasMatrix<Field>& operator()(const Field& F,
								BlasMatrix<Field>& C,
								const BlasMatrix<Field>& A,
								const BlasMatrix<Field>& B) const
		{
			linbox_check( A.rowdim() == B.rowdim());
			linbox_check( C.rowdim() == A.rowdim());
			linbox_check( A.coldim() == B.coldim());
			linbox_check( C.coldim() == A.coldim());
			FFLAS::fsub (F, C.rowdim(), C.coldim(),
				     A.getPointer(), A.getStride(),
				     B.getPointer(), B.getStride(),
				     C.getPointer(), C.getStride());
			return C;
		}
	};

	template<class Field>
	class 	BlasMatrixDomainSubin<Field,BlasMatrix<Field>,BlasMatrix<Field> > {
	public:
		BlasMatrix<Field>& operator()(const Field& F,
								BlasMatrix<Field>& C,
								const BlasMatrix<Field>& B) const
		{
			linbox_check( C.rowdim() == B.rowdim());
			linbox_check( C.coldim() == B.coldim());
			FFLAS::fsubin (F, C.rowdim(), C.coldim(),
				     B.getPointer(), B.getStride(),
				     C.getPointer(), C.getStride());
			return C;
		}
	};

	template<class Field>
	class 	BlasMatrixDomainAddin<Field,BlasMatrix<Field>,BlasMatrix<Field> > {
	public:
		BlasMatrix<Field>& operator()(const Field& F,
								BlasMatrix<Field>& C,
								const BlasMatrix<Field>& B) const
		{
			linbox_check( C.rowdim() == B.rowdim());
			linbox_check( C.coldim() == B.coldim());
			FFLAS::faddin (F, C.rowdim(), C.coldim(),
				     B.getPointer(), B.getStride(),
				     C.getPointer(), C.getStride());
			return C;
		}
	};

	//  general matrix-matrix multiplication and addition with scaling
	// D= beta.C + alpha.A*B
	template<class Field>
	class 	BlasMatrixDomainMulAdd<Field,BlasMatrix<Field>,BlasMatrix<Field>, BlasMatrix<Field> > {
	public:
		BlasMatrix<Field>&
		operator()(const Field                              & F,
			   BlasMatrix<Field>      & D,
			   const typename Field::Element            & beta,
			   const BlasMatrix<Field>& C,
			   const typename Field::Element            & alpha,
			   const BlasMatrix<Field>& A,
			   const BlasMatrix<Field>& B) const
		{
			linbox_check( A.coldim() == B.rowdim());
			linbox_check( C.rowdim() == A.rowdim());
			linbox_check( C.coldim() == B.coldim());
			linbox_check( D.rowdim() == C.rowdim());
			linbox_check( D.coldim() == C.coldim());

			D=C;
			// linbox_check(D.getPointer() != C.getPointer());

			// std::cout << "alpha :" << alpha << std::endl;
			// std::cout << "beta  :" << beta  << std::endl;
			// D.write(std::cout << "Dfgem :=" ) <<','<< std::endl;
			// A.write(std::cout << "Afgem :=" ) <<','<< std::endl;
			// B.write(std::cout << "Bfgem :=" ) <<','<< std::endl;
			FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				      C.rowdim(), C.coldim(), A.coldim(),
				      alpha,
				      A.getPointer(), A.getStride(),
				      B.getPointer(), B.getStride(),
				      beta,
				      D.getPointer(), D.getStride());
			// D.write(std::cout << "Dfgem :=" ) <<','<< std::endl;
			// std::cout << A.getStride() << "," << A.coldim() << std::endl;
			// std::cout << B.getStride() << "," << B.coldim() << std::endl;
			// std::cout << D.getStride() << "," << D.coldim() << std::endl;
			return D;
		}


		BlasMatrix<Field>&
		operator() (const Field                              & F,
			    const typename Field::Element            & beta,
			    BlasMatrix<Field>      & C,
			    const typename Field::Element            & alpha,
			    const BlasMatrix<Field>& A,
			    const BlasMatrix<Field>& B) const
		{
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


	template<class Field>
	class 	BlasMatrixDomainMulAdd<Field,
		BlasMatrix<Field>,
		TransposedBlasMatrix<BlasMatrix<Field> >,
		BlasMatrix<Field> > {
	public:
		BlasMatrix<Field>&
		operator()(const Field                              & F,
			   BlasMatrix<Field>      & D,
			   const typename Field::Element            & beta,
			   const BlasMatrix<Field>& C,
			   const typename Field::Element            & alpha,
			   const TransposedBlasMatrix<BlasMatrix<Field> >& A,
			   const BlasMatrix<Field>& B) const
		{
			linbox_check( A.getMatrix().rowdim() == B.rowdim());
			linbox_check( C.rowdim() == A.getMatrix().coldim());
			linbox_check( C.coldim() == B.coldim());
			linbox_check( D.rowdim() == C.rowdim());
			linbox_check( D.coldim() == C.coldim());

			D=C;

			FFLAS::fgemm( F, FFLAS::FflasTrans, FFLAS::FflasNoTrans,
				      C.rowdim(), C.coldim(), B.rowdim(),
				      alpha,
				      A.getMatrix().getPointer(), A.getMatrix().getStride(),
				      B.getPointer(), B.getStride(),
				      beta,
				      D.getPointer(), D.getStride());

			return D;
		}


		BlasMatrix<Field>&
		operator() (const Field                              & F,
			    const typename Field::Element            & beta,
			    BlasMatrix<Field>      & C,
			    const typename Field::Element            & alpha,
			    const TransposedBlasMatrix<BlasMatrix<Field> >& A,
			    const BlasMatrix<Field>& B) const
		{
			linbox_check( A.getMatrix().rowdim() == B.rowdim());
			linbox_check( C.rowdim() == A.getMatrix().coldim());
			linbox_check( C.coldim() == B.coldim());

			FFLAS::fgemm( F, FFLAS::FflasTrans, FFLAS::FflasNoTrans,
				      C.rowdim(), C.coldim(), B.rowdim(),
				      alpha,
				      A.getMatrix().getPointer(), A.getMatrix().getStride(),
				      B.getPointer(), B.getStride(),
				      beta,
				      C.getPointer(), C.getStride());
			return C;
		}
	};

	template<class Field>
	class 	BlasMatrixDomainMulAdd<Field,
		BlasMatrix<Field>,
		TransposedBlasMatrix<BlasMatrix<Field> >,
		TransposedBlasMatrix<BlasMatrix<Field> > > {
	public:
		BlasMatrix<Field>&
		operator()(const Field                              & F,
			   BlasMatrix<Field>      & D,
			   const typename Field::Element            & beta,
			   const BlasMatrix<Field>& C,
			   const typename Field::Element            & alpha,
			   const TransposedBlasMatrix<BlasMatrix<Field> >& A,
			   const TransposedBlasMatrix<BlasMatrix<Field> >& B) const
		{
			linbox_check( A.getMatrix().rowdim() == B.getMatrix().coldim());
			linbox_check( C.rowdim() == A.getMatrix().coldim());
			linbox_check( C.coldim() == B.getMatrix().rowdim());
			linbox_check( D.rowdim() == C.rowdim());
			linbox_check( D.coldim() == C.coldim());

			D=C;
			// linbox_check(D.getPointer() != C.getPointer());

			FFLAS::fgemm( F, FFLAS::FflasTrans, FFLAS::FflasTrans,
				      C.rowdim(), C.coldim(), A.getMatrix().rowdim(),
				      alpha,
				      A.getMatrix().getPointer(), A.getMatrix().getStride(),
				      B.getMatrix().getPointer(), B.getMatrix().getStride(),
				      beta,
				      D.getPointer(), D.getStride());
			return D;
		}


		BlasMatrix<Field>&
		operator() (const Field                              & F,
			    const typename Field::Element            & beta,
			    BlasMatrix<Field>      & C,
			    const typename Field::Element            & alpha,
			    const TransposedBlasMatrix<BlasMatrix<Field> >& A,
			    const TransposedBlasMatrix<BlasMatrix<Field> >& B) const
		{
			linbox_check( A.getMatrix().rowdim() == B.getMatrix().coldim());
			linbox_check( C.rowdim() == A.getMatrix().coldim());
			linbox_check( C.coldim() == B.getMatrix().rowdim());

			FFLAS::fgemm( F, FFLAS::FflasTrans, FFLAS::FflasTrans,
				      C.rowdim(), C.coldim(), A.getMatrix().rowdim(),
				      alpha,
				      A.getMatrix().getPointer(), A.getMatrix().getStride(),
				      B.getMatrix().getPointer(), B.getMatrix().getStride(),
				      beta,
				      C.getPointer(), C.getStride());
			return C;
		}
	};

	template<class Field>
	class 	BlasMatrixDomainMulAdd<Field,
		BlasMatrix<Field>,
		BlasMatrix<Field>,
		TransposedBlasMatrix<BlasMatrix<Field> > > {
	public:
		BlasMatrix<Field>&
		operator()(const Field                              & F,
			   BlasMatrix<Field>      & D,
			   const typename Field::Element            & beta,
			   const BlasMatrix<Field>& C,
			   const typename Field::Element            & alpha,
			   const BlasMatrix<Field>& A,
			   const TransposedBlasMatrix<BlasMatrix<Field> >& B) const
		{
			linbox_check( A.coldim() == B.getMatrix().coldim());
			linbox_check( C.rowdim() == A.rowdim());
			linbox_check( C.coldim() == B.getMatrix().rowdim());
			linbox_check( D.rowdim() == C.rowdim());
			linbox_check( D.coldim() == C.coldim());

			D=C;
			FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasTrans,
				      C.rowdim(), C.coldim(), A.coldim(),
				      alpha,
				      A.getPointer(), A.getStride(),
				      B.getMatrix().getPointer(), B.getMatrix().getStride(),
				      beta,
				      D.getPointer(), D.getStride());
			return D;
		}


		BlasMatrix<Field>&
		operator() (const Field                              & F,
			    const typename Field::Element            & beta,
			    BlasMatrix<Field>      & C,
			    const typename Field::Element            & alpha,
			    const BlasMatrix<Field>& A,
			    const TransposedBlasMatrix<BlasMatrix<Field> >& B) const
		{
			linbox_check( A.coldim() == B.getMatrix().coldim());
			linbox_check( C.rowdim() == A.rowdim());
			linbox_check( C.coldim() == B.getMatrix().rowdim());

			FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasTrans,
				      C.rowdim(), C.coldim(), A.coldim(),
				      alpha,
				      A.getPointer(), A.getStride(),
				      B.getMatrix().getPointer(), B.getMatrix().getStride(),
				      beta,
				      C.getPointer(), C.getStride());
			return C;
		}
	};

	/*
	 * specialization for Operand1 and Operand3 of type std::vector<Element>
	 * and Operand2 of type BlasMatrix<Field>
	 */

	//  general matrix-vector multiplication and addition with scaling
	// d = beta.c + alpha.A*b
	template<class Field>
	class BlasMatrixDomainMulAdd<Field,std::vector<typename Field::Element>,BlasMatrix<Field>,std::vector<typename Field::Element> > {
	public:
		std::vector<typename Field::Element>& operator() (const Field& F,
								  std::vector<typename Field::Element>& d,
								  const typename Field::Element& beta,
								  const std::vector<typename Field::Element>& c,
								  const typename Field::Element& alpha,
								  const BlasMatrix<Field>& A,
								  const std::vector<typename Field::Element>& b) const
		{
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
								  const BlasMatrix<Field>& A,
								  const std::vector<typename Field::Element>& b) const
		{
			linbox_check( A.coldim() == b.size());
			linbox_check( A.rowdim() == c.size()); //fixed: dpritcha

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

	//  general vector-matrix multiplication and addition with scaling
	// d = beta.c + alpha.a*B -- note order of coldim, rowdim passed to fgemv is switched
	template<class Field>
	class BlasMatrixDomainMulAdd<Field,std::vector<typename Field::Element>,std::vector<typename Field::Element>,BlasMatrix<Field> > {
	public:
		std::vector<typename Field::Element>& operator() (const Field& F,
								  std::vector<typename Field::Element>& d,
								  const typename Field::Element& beta,
								  const std::vector<typename Field::Element>& c,
								  const typename Field::Element& alpha,
								  const std::vector<typename Field::Element>& a,
								  const BlasMatrix<Field>& B) const
		{
			linbox_check( B.rowdim() == a.size());
			linbox_check( B.coldim() == c.size());
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
								  const BlasMatrix<Field>& B) const
		{
			linbox_check( B.rowdim() == a.size());
			linbox_check( B.coldim() == c.size());

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
	 * Specialization for Operand1, Operand2  of type BlasMatrix<Field>
	 * and Operand3 of type BlasPermutation
	 */

	// Matrix permutation product C = A*B
	template<class Field>
	class BlasMatrixDomainMul<Field,BlasMatrix<Field>,BlasMatrix<Field>, BlasPermutation<size_t> > {
	public:
		BlasMatrix<Field>& operator()(const Field& F,
								BlasMatrix<Field>& C,
								const BlasMatrix<Field>& A,
								const BlasPermutation<size_t>& B) const
		{
			C = A;
			return BlasMatrixDomainMulin<Field,BlasMatrix<Field>,BlasPermutation<size_t> >()( F, C, B);
		}
	};

	template<class Field>
	class BlasMatrixDomainMul<Field,BlasMatrix<Field>, BlasPermutation<size_t>,BlasMatrix<Field> > {
	public:
		BlasMatrix<Field>& operator()(const Field& F,
								BlasMatrix<Field>& C,
								const BlasPermutation<size_t>& B,
								const BlasMatrix<Field>& A) const
		{
			C = A;
			return BlasMatrixDomainMulin<Field,BlasMatrix<Field>,BlasPermutation<size_t> >()( F, B, C);
		}
	};

	/*
	 * specialization for Operand1, Operand2  of type BlasMatrix<Field> and Operand3 of type TransposedBlasMatrix<BlasPermutation<size_t> >
	 */

	// Matrix permutation product C = A*B
	template<class Field>
	class BlasMatrixDomainMul<Field,BlasMatrix<Field>,BlasMatrix<Field>, TransposedBlasMatrix<BlasPermutation<size_t> > > {
	public:
		BlasMatrix<Field>& operator()(const Field& F,
								BlasMatrix<Field>& C,
								const BlasMatrix<Field>& A,
								const TransposedBlasMatrix<BlasPermutation<size_t> >& B) const
		{
			C = A;
			return BlasMatrixDomainMulin<Field,BlasMatrix<Field>,TransposedBlasMatrix<BlasPermutation<size_t> > >()( F, C, B);
		}
	};

	template<class Field>
	class BlasMatrixDomainMul<Field,BlasMatrix<Field>, TransposedBlasMatrix<BlasPermutation<size_t> >,BlasMatrix<Field> > {
	public:
		BlasMatrix<Field>& operator()(const Field& F,
								BlasMatrix<Field>& C,
								const TransposedBlasMatrix<BlasPermutation<size_t> >& B,
								const BlasMatrix<Field>& A) const
		{
			C = A;
			return BlasMatrixDomainMulin<Field,BlasMatrix<Field>,TransposedBlasMatrix<BlasPermutation<size_t> > >()( F, B, C);
		}
	};

	/*
	 * specialization for Operand1 of type BlasMatrix<Field> and Operand2 of type BlasPermutation
	 */

	// In-place matrix permutation product
	template<class Field>
	class BlasMatrixDomainMulin<Field,BlasMatrix<Field>, BlasPermutation<size_t> > {
	public:
		BlasMatrix<Field>& operator()( const Field& F,
								 BlasMatrix<Field>& A,
								 const BlasPermutation<size_t>& B) const
		{
			if (B.isIdentity()) return A ;
			linbox_check( A.coldim() >= B.getSize() );
			FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
					A.rowdim(), 0,(int) B.getOrder(),
					A.getPointer(), A.getStride(), B.getPointer() );
			return A;
		}

		BlasMatrix<Field>& operator()( const Field& F,
								 const BlasPermutation<size_t>& B,
								 BlasMatrix<Field>& A) const
		{
			if (B.isIdentity()) return A ;
			linbox_check( A.rowdim() >= B.getSize() );
			FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
					A.coldim(), 0,(int) B.getOrder(), A.getPointer(), A.getStride(), B.getPointer() );
			return A;
		}

	};

	template<class Field>
	class BlasMatrixDomainMulin<Field,BlasMatrix<Field>, TransposedBlasMatrix<BlasPermutation<size_t> > > {
	public:
		BlasMatrix<Field>& operator()( const Field& F,
								 BlasMatrix<Field>& A,
								 const TransposedBlasMatrix<BlasPermutation<size_t> >& B) const
		{
			if (B.getMatrix().isIdentity()) return A ;
			linbox_check( A.coldim() >= B.getMatrix().getSize() );
			FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
					A.rowdim(), 0,(int) B.getMatrix().getOrder(),
					A.getPointer(), A.getStride(), B.getMatrix().getPointer() );
			return A;
		}
		BlasMatrix<Field>& operator()(  const Field& F,
								  const TransposedBlasMatrix<BlasPermutation<size_t> >& B,
								  BlasMatrix<Field>& A) const
		{
			if (B.getMatrix().isIdentity()) return A ;
			linbox_check( A.rowdim() >= B.getMatrix().getSize() );
			FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans,
					A.coldim(), 0,(int) B.getMatrix().getOrder(), A.getPointer(), A.getStride(), B.getMatrix().getPointer() );
			return A;
		}
	};



	/*
	 * specialization for Operand1, Operand2  of type std::vector<Element> and Operand3 of type BlasPermutation
	 */

	// Matrix permutation product C = A*B
	template<class Field>
	class BlasMatrixDomainMul<Field,std::vector< typename Field::Element>,std::vector< typename Field::Element>, BlasPermutation<size_t> > {
	public:
		std::vector< typename Field::Element>& operator()(const Field& F,
								  std::vector< typename Field::Element>& C,
								  const std::vector< typename Field::Element>& A,
								  const BlasPermutation<size_t>& B) const
		{
			C = A;
			return BlasMatrixDomainMulin<Field,std::vector< typename Field::Element>,BlasPermutation<size_t> >()( F, C, B);
		}
	};

	template<class Field>
	class BlasMatrixDomainMul<Field,std::vector< typename Field::Element>, BlasPermutation<size_t>,std::vector< typename Field::Element> > {
	public:
		std::vector< typename Field::Element>& operator()(const Field& F,
								  std::vector< typename Field::Element>& C,
								  const BlasPermutation<size_t>& B,
								  const std::vector< typename Field::Element>& A) const
		{
			C = A;
			return BlasMatrixDomainMulin<Field,std::vector< typename Field::Element>,BlasPermutation<size_t> >()( F, B, C);
		}
	};

	/*
	 * specialization for Operand1, Operand2  of type std::vector<Element> and Operand3 of type TransposedBlasMatrix<BlasPermutation<size_t> >
	 */

	// Matrix permutation product C = A*B
	template<class Field>
	class BlasMatrixDomainMul<Field,std::vector< typename Field::Element>,std::vector< typename Field::Element>, TransposedBlasMatrix<BlasPermutation<size_t> > > {
	public:
		std::vector< typename Field::Element>& operator()(const Field& F,
								  std::vector< typename Field::Element>& C,
								  const std::vector< typename Field::Element>& A,
								  const TransposedBlasMatrix<BlasPermutation<size_t> >& B) const
		{
			C = A;
			return BlasMatrixDomainMulin<Field,std::vector< typename Field::Element>,TransposedBlasMatrix<BlasPermutation<size_t> > >()( F, C, B);
		}
	};

	template<class Field>
	class BlasMatrixDomainMul<Field,std::vector< typename Field::Element>, TransposedBlasMatrix<BlasPermutation<size_t> >,std::vector< typename Field::Element> > {
	public:
		std::vector< typename Field::Element>& operator()(const Field& F,
								  std::vector< typename Field::Element>& C,
								  const TransposedBlasMatrix<BlasPermutation<size_t> >& B,
								  const std::vector< typename Field::Element>& A) const
		{
			C = A;
			return BlasMatrixDomainMulin<Field,std::vector< typename Field::Element>,TransposedBlasMatrix<BlasPermutation<size_t> > >()( F, B, C);
		}
	};

	/*
	 * specialization for Operand1 of type std::vector<Element> and Operand2 of type BlasPermutation
	 */

	// In-place matrix permutation product
	template<class Field>
	class BlasMatrixDomainMulin<Field,std::vector< typename Field::Element>, BlasPermutation<size_t> > {
	public:
		std::vector< typename Field::Element>& operator()( const Field& F,
								   std::vector< typename Field::Element>& A,
								   const BlasPermutation<size_t>& B) const
		{
			if (B.isIdentity()) return A ;
			linbox_check( A.size() == B.getSize() );
			FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
					1, 0,(int) B.getOrder(), &A[0], 1, B.getPointer() );
			return A;
		}

		std::vector< typename Field::Element>& operator()( const Field& F,
								   const BlasPermutation<size_t>& B,
								   std::vector< typename Field::Element>& A) const
		{
			if (B.isIdentity()) return A ;
			linbox_check( A.size() >= B.getSize() );
			FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
					1, 0,(int) B.getOrder(), &A[0], 1, B.getPointer() );
			return A;
		}

	};

	template<class Field>
	class BlasMatrixDomainMulin<Field,std::vector< typename Field::Element>, TransposedBlasMatrix<BlasPermutation<size_t> > > {
	public:
		std::vector< typename Field::Element>& operator()( const Field& F,
								   std::vector< typename Field::Element>& A,
								   const TransposedBlasMatrix<BlasPermutation<size_t> >& B) const
		{
			if (B.getMatrix().isIdentity()) return A ;
			linbox_check( A.size() >= B.getMatrix().getSize() );
			FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
					1, 0,(int) B.getMatrix().getOrder(),
					&A[0], 1, B.getMatrix().getPointer() );
			return A;
		}
		std::vector< typename Field::Element>& operator()(  const Field& F,
								    const TransposedBlasMatrix<BlasPermutation<size_t> >& B,
								    std::vector< typename Field::Element>& A) const
		{
			if (B.getMatrix().isIdentity()) return A ;
			linbox_check( A.size() >= B.getMatrix().getSize() );
			FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans,
					1, 0,(int) B.getMatrix().getOrder(), &A[0], 1, B.getMatrix().getPointer() );
			return A;
		}
	};

	/*
	 * specialization for Operand1 of type BlasMatrix<Field> and Operand2
	 * of type TriangularBlasMatrix<Field>
	 */

	// Matrix/Triangular product C = A*B
	template<class Field>
	class BlasMatrixDomainMul<Field,BlasMatrix<Field>,BlasMatrix<Field>, TriangularBlasMatrix<Field> > {
	public:
		BlasMatrix<Field>& operator()(const Field& F,
								BlasMatrix<Field>& C,
								const BlasMatrix<Field>& A,
								const TriangularBlasMatrix<Field>& B) const
		{
			C = A;
			return BlasMatrixDomainMulin<Field,BlasMatrix<Field>,TriangularBlasMatrix<Field> >()( F, C, B);
		}
	};

	template<class Field>
	class BlasMatrixDomainMul<Field,BlasMatrix<Field>, TriangularBlasMatrix<Field>,BlasMatrix<Field> > {
	public:
		BlasMatrix<Field>& operator()(const Field& F,
								BlasMatrix<Field>& C,
								const TriangularBlasMatrix<Field>& B,
								const BlasMatrix<Field>& A) const
		{
			C = A;
			return BlasMatrixDomainMulin<Field,BlasMatrix<Field>,TriangularBlasMatrix<Field> >()( F, B, C);
		}
	};

	/*
	 * specialization for Operand1 of type BlasMatrix<Field> and Operand2 of type TriangularBlasMatrix<Field>
	 */

	// In-place matrix*triangular matrix product
	template<class Field>
	class BlasMatrixDomainMulin<Field,BlasMatrix<Field>,
	      TriangularBlasMatrix<Field> >{
	public:
		BlasMatrix<Field>& operator()( const Field& F,
								 BlasMatrix<Field>& A,
								 const TriangularBlasMatrix<Field>& B) const
		{
			typename Field::Element one;
			F.init(one, 1UL);
			linbox_check( A.coldim() == B.rowdim() );

			FFLAS::ftrmm( F, FFLAS::FflasRight, (FFLAS::FFLAS_UPLO) (B.getUpLo()),
				      FFLAS::FflasNoTrans,(FFLAS::FFLAS_DIAG) (B.getDiag()),
				      A.rowdim(), A.coldim(), one,
				      B.getPointer(), B.getStride(), A.getPointer(), A.getStride() );
			return A;
		}

		BlasMatrix<Field>& operator()( const Field& F,
								 const TriangularBlasMatrix<Field>& B,
								 BlasMatrix<Field>& A) const
		{
			linbox_check( B.coldim() == A.rowdim() );
			typename Field::Element one;
			F.init(one, 1UL);
			FFLAS::ftrmm( F, FFLAS::FflasLeft, (FFLAS::FFLAS_UPLO)(B.getUpLo()),
				      FFLAS::FflasNoTrans, (FFLAS::FFLAS_DIAG) (B.getDiag()),
				      A.rowdim(), A.coldim(), one,
				      B.getPointer(), B.getStride(),
				      A.getPointer(), A.getStride() );
			return A;
		}
	};


	/*! @internal In-place matrix*triangular matrix product with transpose.
	 */
	template<class Field>
	class BlasMatrixDomainMulin<Field,BlasMatrix<Field>,
	      TransposedBlasMatrix<TriangularBlasMatrix<Field> > >{
	public:
		BlasMatrix<Field>& operator()( const Field& F,
								 BlasMatrix<Field>& A,
								 const TransposedBlasMatrix< TriangularBlasMatrix<Field> >& B) const
		{
			typename Field::Element one;
			F.init(one, 1UL);
			linbox_check( B.getMatrix().coldim() == A.coldim() );

			FFLAS::ftrmm( F, FFLAS::FflasRight,
				      (FFLAS::FFLAS_UPLO)(B.getMatrix().getUpLo()),
				      FFLAS::FflasTrans,
				      (FFLAS::FFLAS_DIAG) (B.getMatrix().getDiag()),
				      A.rowdim(), A.coldim(),
				      one,
				      B.getMatrix().getPointer(), B.getMatrix().getStride(),
				      A.getPointer(), A.getStride() );
			return A;
		}

		BlasMatrix<Field>& operator()( const Field& F,
								 const TransposedBlasMatrix< TriangularBlasMatrix<Field> >& B,
								 BlasMatrix<Field>& A) const
		{
			linbox_check( B.getMatrix().coldim() == A.rowdim() );
			typename Field::Element one;
			F.init(one, 1UL);
			FFLAS::ftrmm( F, FFLAS::FflasLeft,
				      (FFLAS::FFLAS_UPLO) (B.getMatrix().getUpLo()),
				      FFLAS::FflasTrans,
				      (FFLAS::FFLAS_DIAG) (B.getMatrix().getDiag()),
				      A.rowdim(), A.coldim(), one,
				      B.getMatrix().getPointer(), B.getMatrix().getStride(),
				      A.getPointer(), A.getStride() );
			return A;
		}
	};



	/*
	 * specialization for Operand1 of type TriangularBlasMatrix<Field> and Operand2 of type BlasPermutation
	 */

	// Matrix permutation product C = A*B
	template<class Field>
	class BlasMatrixDomainMul<Field,BlasMatrix<Field>,TriangularBlasMatrix<Field>, BlasPermutation<size_t> > {
	public:
		BlasMatrix<Field>& operator()(const Field& F,
								BlasMatrix<Field>& C,
								const TriangularBlasMatrix<Field>& A,
								const BlasPermutation<size_t>& B) const
		{
			C = A;
			return BlasMatrixDomainMulin<Field,BlasMatrix<Field>,BlasPermutation<size_t> >()( F, C, B);
		}
	};

	template<class Field>
	class BlasMatrixDomainMul<Field,BlasMatrix<Field>, BlasPermutation<size_t>,TriangularBlasMatrix<Field> > {
	public:
		BlasMatrix<Field>& operator()(const Field& F,
								BlasMatrix<Field>& C,
								const BlasPermutation<size_t>& B,
								const TriangularBlasMatrix<Field>& A) const
		{
			C = A;
			return BlasMatrixDomainMulin<Field,BlasMatrix<Field>,BlasPermutation<size_t> >()( F, B, C);
		}
	};

	/*
	 * Specialization for Operand of type BlasMatrix<Field>
	 */

	template <class Field>
	class BlasMatrixDomainLeftSolve<Field,BlasMatrix<Field>,BlasMatrix<Field> > {
	public:
		BlasMatrix<Field>& operator() (const Field& F,
								 BlasMatrix<Field>& X,
								 const BlasMatrix<Field>& A,
								 const BlasMatrix<Field>& B) const
		{
			LQUPMatrix<Field> LQUP(A);
			LQUP.left_solve(X,B);
			return X;
		}


		BlasMatrix<Field>& operator() (const Field& F,
								 const BlasMatrix<Field>& A,
								 BlasMatrix<Field>& B) const
		{
			LQUPMatrix<Field> LQUP(A);
			LQUP.left_solve(B);
			return B;
		}
	};

	template <class Field>
	class BlasMatrixDomainRightSolve<Field,BlasMatrix<Field>,BlasMatrix<Field> > {
	public:
		BlasMatrix<Field>& operator() (const Field& F,
								 BlasMatrix<Field>& X,
								 const BlasMatrix<Field>& A,
								 const BlasMatrix<Field>& B) const
		{
			LQUPMatrix<Field> LQUP(A);
			LQUP.right_solve(X,B);
			return X;
		}


		BlasMatrix<Field>& operator() (const Field& F,
								 const BlasMatrix<Field>& A,
								 BlasMatrix<Field>& B) const
		{
			LQUPMatrix<Field> LQUP(A);
			LQUP.right_solve(B);
			return B;
		}

	};

	/*
	 * Specialization for Operand of type std::vector<Element>
	 */

	template <class Field>
	class BlasMatrixDomainLeftSolve<Field, std::vector<typename Field::Element>, BlasMatrix<Field> > {
	public:
		std::vector<typename Field::Element>& operator() (const Field& F,
								  std::vector<typename Field::Element>& X,
								  const BlasMatrix<Field>& A,
								  const std::vector<typename Field::Element>& B) const
		{
			LQUPMatrix<Field> LQUP(A);
			LQUP.left_solve(X,B);
			return X;
		}

		std::vector<typename Field::Element>& operator()(const Field& F,
								 const BlasMatrix<Field>& A,
								 std::vector<typename Field::Element>& B) const
		{
			LQUPMatrix<Field> LQUP(A);
			LQUP.left_solve(B);
			return B;
		}

	};

	template <class Field>
	class BlasMatrixDomainRightSolve<Field, std::vector<typename Field::Element>, BlasMatrix<Field> > {
	public:
		std::vector<typename Field::Element>& operator() (const Field& F,
								  std::vector<typename Field::Element>& X,
								  const BlasMatrix<Field>& A,
								  const std::vector<typename Field::Element>& B) const
		{
			LQUPMatrix<Field> LQUP(A);
			LQUP.right_solve(X,B);
			return X;
		}

		std::vector<typename Field::Element>& operator() (const Field& F,
								  const BlasMatrix<Field>& A,
								  std::vector<typename Field::Element>& B) const
		{
			LQUPMatrix<Field> LQUP(A);
			LQUP.right_solve(B);
			return B;
		}

	};


	/*
	 * ********************************************************
	 * *** Specialization for TriangularBlasMatrix<Field> ***
	 * ********************************************************
	 */


	/*
	 * specialization for Operand of type BlasMatrix<Field>
	 */

	template <class Field>
	class BlasMatrixDomainLeftSolve<Field, BlasMatrix<Field>,TriangularBlasMatrix<Field> > {
	public:
		BlasMatrix<Field>& operator() (const Field& F,
								 BlasMatrix<Field>& X,
								 const TriangularBlasMatrix<Field>& A,
								 const BlasMatrix<Field>& B) const
		{

			linbox_check( X.rowdim() == B.rowdim());
			linbox_check( X.coldim() == B.coldim());

			typename BlasMatrix<Field>::ConstIterator  Biter =   B.Begin();
			typename BlasMatrix<Field>::Iterator       Xiter =   X.Begin();

			for (; Biter != B.End(); ++Biter,++Xiter)
				F.assign(*Xiter,*Biter);

			return (*this)(F,A,X);

		}

		BlasMatrix<Field>& operator() (const Field& F,
								 const TriangularBlasMatrix<Field>& A,
								 BlasMatrix<Field>& B) const
		{
			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.coldim() == B.rowdim());
			typename Field::Element _One;
			F.init(_One,1UL);

			FFLAS::ftrsm( F,
				      FFLAS::FflasLeft, (FFLAS::FFLAS_UPLO) A.getUpLo(),
				      FFLAS::FflasNoTrans,(FFLAS::FFLAS_DIAG) A.getDiag(),
				      A.rowdim(), B.coldim(),
				      _One,A.getPointer(),A.getStride(),
				      B.getPointer(),B.getStride());

			return B;
		}
	};

	template <class Field>
	class BlasMatrixDomainRightSolve<Field,BlasMatrix<Field>, TriangularBlasMatrix<Field> > {
	public:
		BlasMatrix<Field>& operator() (const Field& F,
								 BlasMatrix<Field>& X,
								 const TriangularBlasMatrix<Field>& A,
								 const BlasMatrix<Field>& B) const
		{

			linbox_check( X.rowdim() == B.rowdim());
			linbox_check( X.coldim() == B.coldim());

			typename BlasMatrix<Field>::ConstIterator  Biter =   B.Begin();
			typename BlasMatrix<Field>::Iterator       Xiter =   X.Begin();

			for (; Biter != B.End(); ++Biter,++Xiter)
				F.assign(*Xiter,*Biter);

			return (*this)(F,A,X);
		}

		BlasMatrix<Field>& operator() (const Field& F,
								 const TriangularBlasMatrix<Field>& A,
								 BlasMatrix<Field>& B) const
		{

			linbox_check( A.rowdim() == A.coldim());
			linbox_check( B.coldim() == A.rowdim());
			typename Field::Element _One;
			F.init(_One,1UL);

			FFLAS::ftrsm( F,
				      FFLAS::FflasRight,(FFLAS::FFLAS_UPLO) A.getUpLo(),
				      FFLAS::FflasNoTrans,(FFLAS::FFLAS_DIAG) A.getDiag() ,
				      B.rowdim(), A.coldim(),
				      _One,A.getPointer(),A.getStride(),
				      B.getPointer(),B.getStride());


			return B;
		}
	};

	/*
	 * specialization for Operand of type std::vector<Element>
	 */

	template <class Field>
	class BlasMatrixDomainLeftSolve<Field, std::vector<typename Field::Element>, TriangularBlasMatrix<Field> > {
	public:
		std::vector<typename Field::Element>& operator() (const Field& F,
								  std::vector<typename Field::Element>& x,
								  const TriangularBlasMatrix<Field>& A,
								  const std::vector<typename Field::Element>& b) const
		{

			linbox_check (x.size() == b.size());
			typename std::vector<typename Field::Element>::const_iterator biter = b.begin();
			typename std::vector<typename Field::Element>::iterator       xiter = x.begin();
			for (;biter!=b.end();++biter,++xiter)
				F.assign(*xiter,*biter);

			return (*this)(F,A,x);
		}

		std::vector<typename Field::Element>& operator() (const Field& F,
								  const TriangularBlasMatrix<Field>& A,
								  std::vector<typename Field::Element>& b) const
		{

			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.rowdim() == b.size());

			switch (A.getUpLo()) {
			case LinBoxTag::Upper:
				switch(A.getDiag()) {
				case LinBoxTag::Unit:
					FFLAS::ftrsv( F,
						      FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit,
						      b.size(),A.getPointer(),A.getStride(),&b[0],1);
					break;
				case LinBoxTag::NonUnit:
					FFLAS::ftrsv( F,
						      FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,
						      b.size(),A.getPointer(),A.getStride(),&b[0],1);
					break;
				default:
					throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				}
				break;
			case LinBoxTag::Lower:
				switch(A.getDiag()) {
				case LinBoxTag::Unit:
					FFLAS::ftrsv( F,
						      FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit,
						      b.size(),A.getPointer(),A.getStride(),&b[0],1);
					break;
				case LinBoxTag::NonUnit:
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
	class BlasMatrixDomainRightSolve<Field, std::vector<typename Field::Element>, TriangularBlasMatrix<Field> > {
	public:
		std::vector<typename Field::Element>& operator() (const Field& F,
								  std::vector<typename Field::Element>& x,
								  const TriangularBlasMatrix<Field>& A,
								  const std::vector<typename Field::Element>& b) const
		{

			linbox_check (x.size() == b.size());
			typename std::vector<typename Field::Element>::const_iterator biter = b.begin();
			typename std::vector<typename Field::Element>::iterator       xiter = x.begin();
			for (;biter!=b.end();++biter,++xiter)
				F.assign(*xiter,*biter);

			return (*this)(F,A,x);
		}

		std::vector<typename Field::Element>& operator() (const Field& F,
								  const TriangularBlasMatrix<Field>& A,
								  std::vector<typename Field::Element>& b) const
		{

			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.coldim() == b.size());


			switch (A.getUpLo()) {
			case LinBoxTag::Upper:
				switch(A.getDiag()) {
				case LinBoxTag::Unit:
					FFLAS::ftrsv( F,
						      FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit,
						      b.size(),A.getPointer(),A.getStride(),&b[0],1);
					break;
				case LinBoxTag::NonUnit:
					FFLAS::ftrsv( F,
						      FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit,
						      b.size(),A.getPointer(),A.getStride(),&b[0],1);
					break;
				default:
					throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				}
				break;
			case LinBoxTag::Lower:
				switch(A.getDiag()) {
				case LinBoxTag::Unit:
					FFLAS::ftrsv( F,
						      FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit,
						      b.size(),A.getPointer(),A.getStride(),&b[0],1);
					break;
				case LinBoxTag::NonUnit:
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

	template< class Field, class Polynomial>
	class BlasMatrixDomainMinpoly< Field, Polynomial, BlasMatrix<Field> > {
	public:
		Polynomial& operator() (const Field &F, Polynomial& P, const BlasMatrix<Field>& A) const
		{
			commentator.start ("Modular Dense Minpoly ", "MDMinpoly");

			size_t n = A.coldim();
			linbox_check( n == A.rowdim());
			typename Field::Element * X = new typename Field::Element[n*(n+1)];
			size_t *Perm = new size_t[n];
			for ( size_t i=0; i<n; ++i)
				Perm[i] = 0;
			FFPACK::MinPoly<Field,Polynomial>( F, P, n, A.getPointer(), A.getStride(), X, n, Perm);
			commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "minpoly with " << P.size() << " coefficients" << std::endl;

			delete[] Perm;
			delete[] X;
			commentator.stop ("done",NULL,"MDMinpoly");
			return P;
		}
	};

#if !defined(__INTEL_COMPILER) && !defined(__CUDACC__) && !defined(__clang__)
	template <>
#endif
	template< class Field,  class ContPol >
	class BlasMatrixDomainCharpoly< Field,  ContPol, BlasMatrix<Field> > {
	public:
		ContPol& operator() ( const Field                                	&F,
				      ContPol                     			&P,
				      const BlasMatrix<Field> 	&A) const
		{

			size_t n = A.coldim();
			P.clear();
			linbox_check( n == A.rowdim());
			FFPACK::CharPoly( F, P, n, A.getPointer(), A.getStride());
			return P;
		}
	};

} //end of namespace LinBox

#endif // __LINBOX_blas_matrix_domain_INL

