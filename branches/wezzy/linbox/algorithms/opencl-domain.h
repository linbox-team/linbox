/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/opencl-domain.h
 * Copyright (C) 2011 Matthew Wezowicz, David Saunders
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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/*! @file algorithms/opencl-domain.h
 * @ingroup algorithms
 * @brief NO DOC
 * @warning An <code>OpenCLMatrixDomain<Field></code> should be templated by a
 * Modular<double> or Modular<float> field only.
 */

#ifndef __LINBOX_opencl_matrix_domain_H
#define __LINBOX_opencl_matrix_domain_H

#include <iostream>
#include <vector>
#include <fflas-ffpack/ffpack/ffpack.h>
#include <fflas-ffpack/fflas/fflas.h>
#include "linbox/algorithms/blas-domain.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/matrix/matrix-permutation.h"
#include "linbox/util/debug.h"

#include "CL/cl.hpp"
//#include "helper_functions.hpp" -- For debugging only


namespace LinBox{

	/**
	 *  Interface for all functionnalities provided
	 *  for BlasMatrix using GPUs.
	 *  @internal
	 *  Done through specialization of some of the member funcions
	 *  defined below.  Otherwise, by default the single processor
 	 * BlasMatrixDomain funcions are invoked.
	 */
	template <class Field>
	class OpenCLMatrixDomain {

	public:
		typedef typename Field::Element         Element;

	protected:

		const Field  & _F;
		Element _One;
		Element _Zero;
		Element _MOne;

		//OpenCL specific variables
		cl_context context;
		cl_device_id device;
		cl_command_queue commandQue;
		cl_int errcode;

		//Storage for memory levels
		unsigned long memCapacity;
		unsigned long maxBufferSize;

		//Container type flag
		bool GPUcontainer;
		bool CPUcontainer;
		bool setupCorrect;
		bool doubleSupported;

		//Storage for kernels
		cl_kernel dpKernels[10];
		cl_kernel spKernels[10];

		/**
		 * @internal
		 * Picks the platform used for the container
		 */
		cl_int oclGetPlatformID(cl_platform_id* selected_platform);

		/**
		 * @internal
		 * Picks the device used for the container
		 */
		cl_device_id oclDeviceSelector(cl_int num_devices, cl_device_id* devices);

		/**
		 * @internal
		 * Loads the contents of the specified file into memory
		 * Returns a pointer to a char array and the length of the file
		 */
		char* readFileContents(const char* file_name, int &length);

		/**
		 * @internal
		 * Creates a kernel given a file name and kernel name
		 * Returns the kernel
		 */
		cl_kernel oclCreateKernel(const char* file_name, const char* kernel_name);

		/**
		 * @internal
		 * Initializes the OpenCL compute environment
		 */
		void oclDomainInit();

		/**
		 * @internal
		 * Releases OpenCL cumpute resources
		 */
		void oclDomainTearDown();

		/**
		 * @internal
		 * Checks to see if the memory levels required are possible
		 */
		template<typename T, class Operand1, class Operand2, class Operand3>
		bool oclMemCheck(Operand1 &C, const Operand2 &A, const Operand3 &B) const;

		template<typename T, class Operand1, class Operand2, class Operand3>
		bool oclMemCheck(Operand1& D, const Operand2& A, const Operand3& B, const Operand1& C) const;
		
		template<typename T, class Operand1, class Operand2, class Operand3>
		bool oclMemCheck(Operand1& D, const Operand2& A, const Operand3& B, const Operand1& C, Operand1& Temp) const;

		/**
		 * @internal
		 * OpenCL memory management functions
		 */
		template<class Operand1>
		cl_mem createMatrixBuffer(Operand1 &matrix) const;

		template<class Operand1>
		cl_mem createAndLoadMatrixBuffer(const Operand1 &matrix) const;

		template<class Operand2>
		Operand2& readMatrixBuffer(cl_mem buffer, Operand2 &matrix) const;

		template<typename T, class Operand1>
		cl_mem padMatrix(cl_mem matrixBuffer, int matrixBufferSize,
			int newDimX, const Operand1 &matrix) const;

		template<typename T, class Operand1>
		Operand1& depadMatrix(cl_mem matrixBuffer, int matrixBufferSize,
			int outputSize, int newDimX, Operand1& matrix) const;

	public:

		//! Constructor of OpenCLDomain.

		OpenCLMatrixDomain (const Field& F ) :
			_F(F)
		{
			F.init(_One,1UL);
			F.init(_Zero,0UL);
			F.init(_MOne,-1);
#ifndef NDEBUG
			if (!Givaro::probab_prime(F.characteristic())) {
				std::cout << " *** WARNING *** "                                           << std::endl;
				std::cout << " You are using a BLAS Domain where your field is not prime " << std::endl;
			}
#endif

			//Initialize OpenCL environment
			oclDomainInit(); //TODO -- make configure time instead of run time

		}

		//! Copy constructor
		OpenCLMatrixDomain (const OpenCLMatrixDomain<Field> & BMD) :
			_F(BMD._F), _One(BMD._One), _Zero(BMD._Zero), _MOne(BMD._MOne)
		{
#ifndef NDEBUG
			if (!Givaro::probab_prime(_F.characteristic())) {
				std::cout << " *** WARNING *** "                                           << std::endl;
				std::cout << " You are using a BLAS Domain where your field is not prime " << std::endl;
			}
#endif

			//Initialize OpenCL environment
			oclDomainInit(); //TODO -- make configure time instead of run time

		}

		//! Deconstructor
		~OpenCLMatrixDomain(){
			oclDomainTearDown();
		}


		//! Field accessor
		const Field& field() const
		{
			return _F;
		}


		/*
		 * Basics operation available matrix respecting BlasMatrix interface
		 */

		//! multiplication.
		//! C = A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& mul(Operand1& C, const Operand2& A, const Operand3& B) const
		{
			return BlasMatrixDomainMul<Field,Operand1,Operand2,Operand3>()(_F,C,A,B);
		}

		//! addition.
		//! C = A+B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& add(Operand1& C, const Operand2& A, const Operand3& B) const
		{
			return BlasMatrixDomainAdd<Field,Operand1,Operand2,Operand3>()(_F,C,A,B);
		}

		//! copy.
		//! B = A
		template <class Operand1, class Operand2>
		Operand1& copy(Operand1& B, const Operand2& A) const
		{
			return BlasMatrixDomainCopy<Field,Operand1,Operand2>()(_F,B,A);
		}

		//! substraction
		//! C = A-B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& sub(Operand1& C, const Operand2& A, const Operand3& B) const
		{
			return BlasMatrixDomainSub<Field,Operand1,Operand2,Operand3>()(_F,C,A,B);
		}

		//! substraction (in place)
		//! C -= B
		template <class Operand1, class Operand3>
		Operand1& subin(Operand1& C, const Operand3& B) const
		{
			return BlasMatrixDomainSubin<Field,Operand1,Operand3>()(_F,C,B);
		}

		//! addition (in place)
		//! C += B
		template <class Operand1, class Operand3>
		Operand1& addin(Operand1& C, const Operand3& B) const
		{
			return BlasMatrixDomainAddin<Field,Operand1,Operand3>()(_F,C,B);
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
			return BlasMatrixDomainMulin<Field,Operand1,Operand2>()(_F,A,B);
		}

		//! In place multiplication.
		//! B = A*B
		template <class Operand1, class Operand2>
		Operand2& mulin_right(const Operand1& A, Operand2& B ) const
		{
			return BlasMatrixDomainMulin<Field,Operand2,Operand1>()(_F,A,B);
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
			return BlasMatrixDomainMulAdd<Field,Operand1,Operand2,Operand3>()(_F,D,beta,C,alpha,A,B);
		}

		//! muladdin.
		//! C= beta.C + alpha.A*B.
		template <class Operand1, class Operand2, class Operand3>
		Operand1& muladdin(const Element& beta, Operand1& C,
				   const Element& alpha, const Operand2& A, const Operand3& B) const
		{
			return BlasMatrixDomainMulAdd<Field,Operand1,Operand2,Operand3>()(_F,beta,C,alpha,A,B);
		}


		/*!
		 * @name Solutions available for matrix respecting BlasMatrix interface
		 */
		//@{

		//! Inversion
		template <class Matrix>
		Matrix& inv( Matrix &Ainv, const Matrix &A) const
		{
			BlasMatrixDomainInv<Field,Matrix>()(_F,Ainv,A);
			return Ainv;
		}

		//! Inversion (in place)
		template <class Matrix>
		Matrix& invin( Matrix &Ainv, Matrix &A) const
		{
			BlasMatrixDomainInv<Field,Matrix>()(_F,Ainv,A);
			return Ainv;
		}

		//! Inversion (the matrix A is modified)
		template <class Matrix>
		Matrix& invin(Matrix &A) const
		{
			Matrix tmp(A.rowdim(), A.coldim());
			tmp = A;
			BlasMatrixDomainInv<Field,Matrix>()(_F,A,tmp);
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
			nullity = BlasMatrixDomainInv<Field,Matrix>()(_F,Ainv,A);
			return Ainv;
		}

		//! Inversion (the matrix A is modified) w singular check
		template <class Matrix>
		Matrix& invin( Matrix &Ainv, Matrix &A, int& nullity) const
		{
			nullity = BlasMatrixDomainInv<Field,Matrix>()(_F,Ainv,A);
			return Ainv;
		}

		//! Rank
		template <class Matrix>
		unsigned int rank(const Matrix &A) const
		{
			return BlasMatrixDomainRank<Field,Matrix>()(_F,A);
		}

		//! in-place Rank (the matrix is modified)
		template <class Matrix>
		unsigned int rankin(Matrix &A) const
		{
			return BlasMatrixDomainRank<Field, Matrix>()(_F,A);
		}

		//! determinant
		template <class Matrix>
		Element det(const Matrix &A) const
		{
			return BlasMatrixDomainDet<Field, Matrix>()(_F,A);
		}

		//! in-place Determinant (the matrix is modified)
		template <class Matrix>
		Element detin(Matrix &A) const
		{
			return BlasMatrixDomainDet<Field, Matrix>()(_F,A);
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
			return BlasMatrixDomainLeftSolve<Field,Operand,Matrix>()(_F,X,A,B);
		}

		//! linear solve with matrix right hand side, the result is stored in-place in B.
		//! @pre A must be square
		//! AX=B , (B<-X)
		template <class Operand,class Matrix>
		Operand& left_solve (const Matrix& A, Operand& B) const
		{
			return BlasMatrixDomainLeftSolve<Field,Operand,Matrix>()(_F,A,B);
		}

		//! linear solve with matrix right hand side.
		//! XA=B
		template <class Operand, class Matrix>
		Operand& right_solve (Operand& X, const Matrix& A, const Operand& B) const
		{
			return BlasMatrixDomainRightSolve<Field,Operand,Matrix>()(_F,X,A,B);
		}

		//! linear solve with matrix right hand side, the result is stored in-place in B.
		//! @pre A must be square
		//! XA=B , (B<-X)
		template <class Operand, class Matrix>
		Operand& right_solve (const Matrix& A, Operand& B) const
		{
			return BlasMatrixDomainRightSolve<Field,Operand,Matrix>()(_F,A,B);
		}

		//! minimal polynomial computation.
		template <class Polynomial, class Matrix>
		Polynomial& minpoly( Polynomial& P, const Matrix& A ) const
		{
			return BlasMatrixDomainMinpoly<Field, Polynomial, Matrix>()(_F,P,A);
		}

		//! characteristic polynomial computation.
		template <class Polynomial,  class Matrix >
		Polynomial& charpoly( Polynomial& P, const Matrix& A ) const
		{

			commentator.start ("Modular Dense Charpoly ", "MDCharpoly");
			std::list<Polynomial> P_list;
			P_list.clear();
			BlasMatrixDomainCharpoly<Field, std::list<Polynomial>, Matrix >()(_F,P_list,A);


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
			commentator.stop ("done", NULL, "MDCharpoly");

			return P;
		}

		//! characteristic polynomial computation.
		template <class Polynomial, class Matrix >
		std::list<Polynomial>& charpoly( std::list<Polynomial>& P, const Matrix& A ) const
		{
			return BlasMatrixDomainCharpoly<Field, std::list<Polynomial>, Matrix >()(_F,P,A);
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
				_F.assign(res[i],_Zero);
			for ( i=0;i<P1.size();i++)
				for ( j=0;j<P2.size();j++)
					_F.axpyin(res[i+j],P1[i],P2[j]);
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
					if (!_F.areEqual(A.getEntry(i,j),B.getEntry(i,j))) //!@bug use refs
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
					if (!_F.isZero(A.getEntry(i,j))) //!@bug use refs
						return false ;
			return true ;
		}

		template<class Matrix1>
		bool isIdentity(const Matrix1 & A)
		{
			if (A.rowdim() != A.coldim())
				return false ;
			for (size_t i = 0 ; i < A.rowdim() ; ++i)
				if (!_F.isOne(A.getEntry(i,i)))
					return false;

			for (size_t i = 0 ; i < A.rowdim() ; ++i)
				for (size_t j = 0 ; j < i ; ++j)
					if (!_F.isZero(A.getEntry(i,j))) //!@bug use refs
						return false ;
			for (size_t i = 0 ; i < A.rowdim() ; ++i)
				for (size_t j = i+1 ; j < A.coldim() ; ++j)
					if (!_F.isZero(A.getEntry(i,j))) //!@bug use refs
						return false ;
			return true ;
		}

		template<class Matrix1>
		bool isIdentityGeneralized(const Matrix1 & A)
		{
			size_t mn = std::min(A.rowdim(),A.coldim());
			for (size_t i = 0 ; i < mn ; ++i)
				if (!_F.isOne(A.getEntry(i,i)))
					return false;

			for (size_t i = 0 ; i < A.rowdim() ; ++i)
				for (size_t j = 0 ; j < std::min(i,mn) ; ++j)
					if (!_F.isZero(A.getEntry(i,j))) //!@bug use refs
						return false ;
			for (size_t i = 0 ; i < A.rowdim() ; ++i)
				for (size_t j = i+1 ; j < A.coldim() ; ++j)
					if (!_F.isZero(A.getEntry(i,j))) //!@bug use refs
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
			return A.write (os, _F);
		}

		template <class Matrix>
		inline std::ostream &write (std::ostream &os, const Matrix &A, bool maple_format) const
		{
			return A.write (os, _F, maple_format);
		}

		/** Read matrix
		 * @param  is  Input stream from which matrix is read.
		 * @param  A   Matrix.
		 * @returns reference to is.
		 */
		template <class Matrix>
		inline std::istream &read (std::istream &is, Matrix &A) const
		{
			return A.read (is, _F);
		}

	}; /* end of class OpenCLMatrixDomain */

} /* end of namespace LinBox */

#include "linbox/algorithms/opencl-domain-setup.inl"
#include "linbox/algorithms/opencl-domain-memory.inl"
#include "linbox/algorithms/opencl-domain.inl"

#endif /* __LINBOX_opencl_matrix_domain_H */

