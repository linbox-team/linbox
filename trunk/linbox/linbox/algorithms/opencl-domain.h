/* linbox/algorithms/opencl-domain.h
 * Copyright (C) 2011      David Saunders
 *               2011-2012 Matthew Wezowicz
 *
 * Written by Matthew Wezowicz <mwezz@udel.edu>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file algorithms/opencl-domain.h
 * @ingroup algorithms
 * @brief NO DOC
 * @warning An <code>OpenCLMatrixDomain<Field></code> should be templated by a
 * Modular<double> or Modular<float> field only.
 */

#ifndef __LINBOX_opencl_matrix_domain_H
#define __LINBOX_opencl_matrix_domain_H

#include <vector>
#include <iostream>
#include <pthread.h>

#include "linbox/algorithms/blas-domain.h"
#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#ifdef __LINBOX_HAVE_OCL
#include "CL/cl.h"
#endif

namespace LinBox{

	/**
	 * Generic submatrix view adapter used internally in the OpenCLMatrixDomain
	 */
	template <class _Matrix>
	class SubmatrixAdapter{
	public:
		//Access to underlying types
		typedef typename _Matrix::Field     Field;
		typedef typename Field::Element     Element;
		typedef SubmatrixAdapter<_Matrix>   Self_t;

	private:
		_Matrix* _Mat;  //!< Parent Matrix (ie raw vector)
		size_t _row;    //!< row dimension of Submatrix
		size_t _col;    //!< col dimension of Submatrix
		size_t _r0;     //!< upper left corner row of Submatrix in \p _Mat
		size_t _c0;     //!< upper left corner row of Submatrix in \p _Mat
		size_t _stride; //!< number of columns in \p _Mat (or stride of \p _Mat)
		size_t _off;    //!< offset from start of parent matrix

	public:
		/** NULL constructor.  */
		SubmatrixAdapter() :
			_Mat(NULL),
			_row(0),
			_col(0),
			_r0(0),
			_c0(0),
			_stride(0),
			_off(0){}

		/** Constructor from an existing @refMatrix
		 * \param M Pointer to @ref Matrix of which to construct submatrix
		 */
		SubmatrixAdapter(const _Matrix& M) :
			_Mat(&(const_cast<_Matrix&>(M))),
			_row(M.rowdim()),
			_col(M.coldim()),
			_r0(0),
			_c0(0),
			_stride(M.coldim()),
			_off(0){}

		/** Constructor from an existing @ref Matrix and dimensions.
		 * \param M Pointer to @ref Matrix of which to construct submatrix
		 * \param row Starting row
		 * \param col Starting column
		 * \param Rowdim Row dimension
		 * \param Coldim Column dimension
		 */
		SubmatrixAdapter(
			const _Matrix& M,
			size_t row,
			size_t col,
			size_t Rowdim,
			size_t Coldim) :
			//Init list starts here
			_Mat(&(const_cast<_Matrix&>(M))),
			_row(Rowdim),
			_col(Coldim),
			_r0(row),
			_c0(col),
			_stride(M.coldim()),
			_off(row * _stride + col){}

		/** Constructor from an existing @ref SubmatrixAdapter
		 * \param SM Pointer to @ref SubmatrixAdapter of which to construct submatrix
		 */
		SubmatrixAdapter(const SubmatrixAdapter<_Matrix>& SM) :
			_Mat(SM._Mat),
			_row(SM._row),
			_col(SM._col),
			_r0(SM._r0),
			_c0(SM._c0),
			_stride(SM._stride),
			_off(SM._off){}

		/** Constructor from an existing submatrix and dimensions
		 * @param SM Constant reference to SubmatrixAdapter from which to
		 *           construct submatrix
		 * @param rowbeg Starting row
		 * @param colbeg Starting column
		 * @param Rowdim Row dimension
		 * @param Coldim Column dimension
		 */
		SubmatrixAdapter(
			const SubmatrixAdapter<_Matrix>& SM,
			size_t row,
			size_t col,
			size_t Rowdim,
			size_t Coldim) :
			//Init list starts here
			_Mat(SM._Mat),
			_row(Rowdim),
			_col(Coldim),
			_r0(SM._r0 + row),
			_c0(SM._c0 + col),
			_stride(SM._stride),
			_off(SM._off + (row * _stride + col)){}

		/** Get the number of rows in the matrix
		 * @return Number of rows in matrix
		 */
		size_t rowdim() const{
			return _row;
		}

		/** Get the number of columns in the matrix
		 * @return Number of columns in matrix
		 */
		size_t coldim() const{
			return _col;
		}

		/*! Get the stride of the matrix.
		 * @return stride of submatrix (number of cols of parent matrix)
		 */
		size_t getStride() const{
			return _stride;
		}

		/** Set the entry at (i, j).
		 * @param i Row index of entry, 0...rowdim() - 1
		 * @param j Column index of entry, 0...coldim() - 1
		 * @param a_ij Element to set
		 */
		void setEntry(size_t i, size_t j, const Element& a_ij){
			_Mat->setEntry(_r0 + i, _c0 + j, a_ij);
		}

		/** Get a writeable reference to an entry in the matrix.
		 * @param i Row index of entry,  0...rowdim() - 1
		 * @param j Column index of entry, 0...coldim() - 1
		 * @return Reference to matrix entry
		 */
		Element& refEntry(size_t i, size_t j){
			_Mat->refEntry(_r0 + i, _c0 + j);
		}

		/** Get a read-only individual entry from the matrix.
		 * @param i Row index of entry,  0...rowdim() - 1
		 * @param j Column index of entry, 0...coldim() - 1
		 * @return Const reference to matrix entry
		 */
		const Element& getEntry(size_t i, size_t j) const{
			return _Mat->getEntry(_r0 + i, _c0 + j);
		}

		/** Get an entry and store it in the given value.
		 * This form is more in the Linbox style and is provided for interface
		 * compatibility with other parts of the library
		 * @param x Element in which to store result
		 * @param i Row index of entry,  0...rowdim() - 1
		 * @param j Column index of entry, 0...coldim() - 1
		 * @return Reference to x
		 */
		Element& getEntry(Element& x, size_t i, size_t j){
			return _Mat->getEntry(x, _r0 + i, _c0 + j);
		}

		/** Access the parent matrix
		 * @return Reference to _Mat
		 */
		_Matrix& getMatrix(){
			return *_Mat;
		}
	};

	/**
	 * Interface for all functionnalities provided
	 * for BlasMatrix using GPUs.
	 * @internal
	 * Done through specialization of some of the member funcions
	 * defined below.  Otherwise, by default the single processor
 	 * BlasMatrixDomain funcions are invoked.
	 */
	template <class Field_>
	class OpenCLMatrixDomain {

	public:
		typedef Field_                          Field;
		typedef typename Field::Element         Element;
		typedef BlasMatrix<Field>               Matrix;
#ifdef __LINBOX_HAVE_OCL
		friend class OpenCLMatrixDomainFactory;
#endif

	protected:

		const Field& _F;

#ifdef __LINBOX_HAVE_OCL
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

		//Storage for kernels and flags for availability
		cl_kernel dpKernels[20];
		bool dpKernelsAvailable[20];
		cl_kernel spKernels[20];
		bool spKernelsAvailable[20];

		//ID number assigned by OpenCLMatrixDomainFactory
		unsigned int IDnum;

		//Mutex
		pthread_mutex_t* deviceLock;

		/**
		 * @internal
		 * Initializes the OpenCL compute environment
		 */
		void oclMatrixDomainAcquire();

		/**
		 * @internal
		 * Releases OpenCL cumpute resources
		 */
		void oclMatrixDomainRelease(unsigned int IDnum);

		/**
		 * @internal
		 * Checks to see if the memory levels required are possible
		 */
		template <class Operand1, class Operand2, class Operand3>
		bool oclMemCheck(
			Operand1& D,
			const Operand2& A,
			const Operand3& B,
			const Operand1& C) const;

		template <class Operand1>
		bool oclMemCheck(
			Operand1& D,
			Operand1& A,
			Operand1& B,
			Operand1& C) const;

		/**
		 * @internal
		 * OpenCL memory management functions
		 */
		template <typename T, class Operand1>
		cl_mem oclCreateMatrixBuffer(Operand1 &matrix) const;

		template <typename T, class Operand1>
		cl_mem oclCreateAndLoadMatrixBuffer(const Operand1 &matrix) const;

		template <typename T, class Operand2>
		Operand2& oclReadMatrixBuffer(cl_mem buffer, Operand2 &matrix) const;

		template <typename T, class Operand1>
		cl_mem oclPadMatrix(
			cl_mem matrixBuffer,
			int matrixBufferSize,
			int newDimX,
			const Operand1 &matrix) const;

		template <typename T, class Operand1>
		Operand1& oclDepadMatrix(
			cl_mem matrixBuffer,
			int matrixBufferSize,
			int outputSize,
			int newDimX,
			Operand1& matrix) const;

		/**
		 * @internal
		 * Functions to call the passed kernel on the passed buffers
		 */
		template <typename T, typename U>
		void oclCallKernel(
			cl_mem bufferC,
			cl_mem bufferA,
			cl_mem bufferB,
			int widthA,
			int heightA,
			int widthB,
			T p,
			cl_kernel selectedKernel) const;

		template <typename T, typename U>
		void oclCallKernel(
			cl_mem bufferD,
			cl_mem bufferA,
			cl_mem bufferB,
			cl_mem bufferC,
			int widthA,
			int heightA,
			int widthB,
			T p,
			cl_kernel selectedKernel) const;

		template <typename T, typename U>
		void oclCallKernel(
			cl_mem bufferD,
			cl_mem bufferA,
			cl_mem bufferB,
			cl_mem bufferC,
			T alpha,
			T beta,
			int widthA,
			int heightA,
			int widthB,
			T p,
			cl_kernel selectedKernel) const;

		/**
		 * @internal
		 * Functions to partition the matrices into submatrix views
		 */
		template <class Operand1, class Operand2, class Operand3>
		std::vector<int> oclPartition(
			Operand1& C,
			const Operand2& A,
			const Operand3& B,
			std::vector<SubmatrixAdapter<Operand1> >& VC,
			std::vector<SubmatrixAdapter<Operand2> >& VA,
			std::vector<SubmatrixAdapter<Operand3> >& VB) const;

		template <class Operand1, class Operand2, class Operand3>
		std::vector<int> oclPartition(
			Operand1& D,
			const Operand2& A,
			const Operand3& B,
			const Operand1& C,
			std::vector<SubmatrixAdapter<Operand1> >& VD,
			std::vector<SubmatrixAdapter<Operand2> >& VA,
			std::vector<SubmatrixAdapter<Operand3> >& VB,
			std::vector<SubmatrixAdapter<Operand1> >& VC) const;

		void printClErrstring(cl_int err) const;
#else
		bool setupCorrect;
#endif
	public:

		//! Constructor of OpenCLDomain.
		OpenCLMatrixDomain(const Field& F ) : _F(F), setupCorrect(false){

#ifndef NDEBUG
			if(!Givaro::probab_prime(_F.characteristic())){
				std::cout << " *** WARNING *** " << std::endl;
				std::cout << " You are using a OpenCL Matrix Domain"
				          << " where your field is not prime "
				          << std::endl;
			}
#endif

#ifdef __LINBOX_HAVE_OCL
			//Initialize OpenCL environment
			oclMatrixDomainAcquire();
#endif
		}

		//! Copy constructor
		OpenCLMatrixDomain(const OpenCLMatrixDomain<Field> & OMD) :
			_F(OMD._F),
			setupCorrect(false){

#ifndef NDEBUG
			if(!Givaro::probab_prime(_F.characteristic())){
				std::cout << " *** WARNING *** " << std::endl;
				std::cout << " You are using a OpenCL Matrix Domain"
				          << " where your field is not prime "
				          << std::endl;
			}
#endif

#ifdef __LINBOX_HAVE_OCL
			//Initialize OpenCL environment
			oclMatrixDomainAcquire();
#endif
		}

		//! Deconstructor
		~OpenCLMatrixDomain(){
#ifdef __LINBOX_HAVE_OCL
			oclMatrixDomainRelease(IDnum);
#endif
		}

		//! Field accessor
		const Field& field() const{
			return _F;
		}

		/*
		 * Basics operation available matrix respecting BlasMatrix interface
		 */

		//! multiplication.
		//! C = A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& mul(Operand1& C, const Operand2& A, const Operand3& B) const{
			return BlasMatrixDomainMul<Field,Operand1,Operand2,Operand3>()(_F,C,A,B);
		}

		//! addition.
		//! C = A+B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& add(Operand1& C, const Operand2& A, const Operand3& B) const{
			return BlasMatrixDomainAdd<Field,Operand1,Operand2,Operand3>()(_F,C,A,B);
		}

		//! copy.
		//! B = A
		template <class Operand1, class Operand2>
		Operand1& copy(Operand1& B, const Operand2& A) const{
			return BlasMatrixDomainCopy<Field,Operand1,Operand2>()(_F,B,A);
		}

		//! substraction
		//! C = A-B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& sub(Operand1& C, const Operand2& A, const Operand3& B) const{
			return BlasMatrixDomainSub<Field,Operand1,Operand2,Operand3>()(_F,C,A,B);
		}

		//! substraction (in place)
		//! C -= B
		template <class Operand1, class Operand3>
		Operand1& subin(Operand1& C, const Operand3& B) const{
			return BlasMatrixDomainSubin<Field,Operand1,Operand3>()(_F,C,B);
		}

		//! addition (in place)
		//! C += B
		template <class Operand1, class Operand3>
		Operand1& addin(Operand1& C, const Operand3& B) const{
			return BlasMatrixDomainAddin<Field,Operand1,Operand3>()(_F,C,B);
		}

		//! multiplication with scaling.
		//! C = alpha.A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& mul(
			Operand1& C,
			const Element& alpha,
			const Operand2& A,
			const Operand3& B) const{

			return muladdin(_F.zero,C,alpha,A,B);
		}

		//! In place multiplication.
		//! A = A*B
		template <class Operand1, class Operand2>
		Operand1& mulin_left(Operand1& A, const Operand2& B) const{
			return BlasMatrixDomainMulin<Field,Operand1,Operand2>()(_F,A,B);
		}

		//! In place multiplication.
		//! B = A*B
		template <class Operand1, class Operand2>
		Operand2& mulin_right(const Operand1& A, Operand2& B) const{
			return BlasMatrixDomainMulin<Field,Operand2,Operand1>()(_F,A,B);
		}

		//! axpy.
		//! D = A*B + C
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axpy(
			Operand1& D,
			const Operand2& A,
			const Operand3& B,
			const Operand1& C) const{

			return muladd(D,_F.one,C,_F.one,A,B);
		}

		//! axpyin.
		//! C += A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axpyin(Operand1& C, const Operand2& A, const Operand3& B) const{
			return muladdin(_F.one,C,_F.one,A,B);
		}

		//! maxpy.
		//! D = C - A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& maxpy(
			Operand1& D,
			const Operand2& A,
			const Operand3& B,
			const Operand1& C) const{

			return muladd(D,_F.one,C,_F.mOne,A,B);
		}

		//! maxpyin.
		//! C -= A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& maxpyin(Operand1& C, const Operand2& A, const Operand3& B) const{
			return muladdin(_F.one,C,_F.mOne,A,B);
		}

		//! axmy.
		//! D= A*B - C
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axmy(
			Operand1& D,
			const Operand2& A,
			const Operand3& B,
			const Operand1& C) const{

			return muladd(D,_F.mOne,C,_F.one,A,B);
		}

		//! axmyin.
		//! C = A*B - C
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axmyin(Operand1& C, const Operand2& A, const Operand3& B) const{
			return muladdin(_F.mOne,C,_F.one,A,B);
		}

		//!  general matrix-matrix multiplication and addition with scaling.
		//! D= beta.C + alpha.A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& muladd(
			Operand1& D,
			const Element& beta,
			const Operand1& C,
			const Element& alpha,
			const Operand2& A,
			const Operand3& B) const{

			return BlasMatrixDomainMulAdd<Operand1,Operand2,Operand3>()(
				_F,
				D,
				beta,
				C,
				alpha,
				A,
				B);
		}

		//! muladdin.
		//! C= beta.C + alpha.A*B.
		template <class Operand1, class Operand2, class Operand3>
		Operand1& muladdin(
			const Element& beta,
			Operand1& C,
			const Element& alpha,
			const Operand2& A,
			const Operand3& B) const{

			return BlasMatrixDomainMulAdd<Operand1,Operand2,Operand3>()(
				_F,
				beta,
				C,
				alpha,
				A,
				B);
		}

		/*!
		 * @name Solutions available for matrix respecting BlasMatrix interface
		 */
		//@{

		//! Inversion
		template <class Matrix>
		Matrix& inv( Matrix &Ainv, const Matrix &A) const{
			BlasMatrixDomainInv<Field,Matrix,Matrix>()(_F,Ainv,A);
			return Ainv;
		}

		//! Inversion (in place)
		template <class Matrix>
		Matrix& invin( Matrix &Ainv, Matrix &A) const{
			BlasMatrixDomainInv<Field,Matrix,Matrix>()(_F,Ainv,A);
			return Ainv;
		}

		//! Inversion (the matrix A is modified)
		template <class Matrix>
		Matrix& invin(Matrix &A) const{
			Matrix tmp(A.rowdim(), A.coldim());
			tmp = A;
			BlasMatrixDomainInv<Field,Matrix,Matrix>()(_F,A,tmp);
			return A;
		}

		/*! Division.
		 * C = A B^{-1}  ==>  C . B = A
		 */
		template <class Matrix>
		Matrix& div(Matrix &C, const Matrix &A, const Matrix &B) const{
			return this->right_solve(C,B,A);
		}


		//! Inversion w singular check
		template <class Matrix>
		Matrix& inv(Matrix &Ainv, const Matrix &A, int& nullity) const{
			nullity = BlasMatrixDomainInv<Field,Matrix,Matrix>()(_F,Ainv,A);
			return Ainv;
		}

		//! Inversion (the matrix A is modified) w singular check
		template <class Matrix>
		Matrix& invin(Matrix &Ainv, Matrix &A, int& nullity) const{
			nullity = BlasMatrixDomainInv<Field,Matrix,Matrix>()(_F,Ainv,A);
			return Ainv;
		}

		//! Rank
		template <class Matrix>
		unsigned int rank(const Matrix &A) const{
			return BlasMatrixDomainRank<Field,Matrix>()(_F,A);
		}

		//! in-place Rank (the matrix is modified)
		template <class Matrix>
		unsigned int rankin(Matrix &A) const{
			return BlasMatrixDomainRank<Field, Matrix>()(_F,A);
		}

		//! determinant
		template <class Matrix>
		Element det(const Matrix &A) const{
			return BlasMatrixDomainDet<Field, Matrix>()(_F,A);
		}

		//! in-place Determinant (the matrix is modified)
		template <class Matrix>
		Element detin(Matrix &A) const{
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
		Operand& left_solve (Operand& X, const Matrix& A, const Operand& B) const{
			return BlasMatrixDomainLeftSolve<Field,Operand,Matrix>()(_F,X,A,B);
		}

		//! linear solve with matrix right hand side, the result is stored in-place in B.
		//! @pre A must be square
		//! AX=B , (B<-X)
		template <class Operand,class Matrix>
		Operand& left_solve (const Matrix& A, Operand& B) const{
			return BlasMatrixDomainLeftSolve<Field,Operand,Matrix>()(_F,A,B);
		}

		//! linear solve with matrix right hand side.
		//! XA=B
		template <class Operand, class Matrix>
		Operand& right_solve (Operand& X, const Matrix& A, const Operand& B) const{
			return BlasMatrixDomainRightSolve<Field,Operand,Matrix>()(_F,X,A,B);
		}

		//! linear solve with matrix right hand side, the result is stored in-place in B.
		//! @pre A must be square
		//! XA=B , (B<-X)
		template <class Operand, class Matrix>
		Operand& right_solve (const Matrix& A, Operand& B) const{
			return BlasMatrixDomainRightSolve<Field,Operand,Matrix>()(_F,A,B);
		}

		//! minimal polynomial computation.
		template <class Polynomial, class Matrix>
		Polynomial& minpoly( Polynomial& P, const Matrix& A ) const{
			return BlasMatrixDomainMinpoly<Field, Polynomial, Matrix>()(_F,P,A);
		}

		//! characteristic polynomial computation.
		template <class Polynomial,  class Matrix >
		Polynomial& charpoly( Polynomial& P, const Matrix& A ) const{

			commentator().start ("Modular Dense Charpoly ", "MDCharpoly");
			std::list<Polynomial> P_list;
			P_list.clear();
			BlasMatrixDomainCharpoly<Field, std::list<Polynomial>, Matrix >()(
				_F,
				P_list,
				A);

			Polynomial tmp(A.rowdim() + 1);
			typename std::list<Polynomial>::iterator it = P_list.begin();
			P = *(it++);
			while(it != P_list.end()){
				// Waiting for an implementation of a domain of polynomials
				mulpoly(tmp, P, *it);
				P = tmp;
				//delete &(*it);
				++it;
			}
			commentator().stop ("done", NULL, "MDCharpoly");

			return P;
		}

		//! characteristic polynomial computation.
		template <class Polynomial, class Matrix >
		std::list<Polynomial>& charpoly(
			std::list<Polynomial>& P,
			const Matrix& A ) const{

			return BlasMatrixDomainCharpoly<
				Field,
				std::list<Polynomial>,
				Matrix >()(_F,P,A);
		}


		//private:
		//! @todo Temporary: waiting for an implementation of a domain of polynomial
		template<class Polynomial>
		Polynomial& mulpoly(
			Polynomial &res,
			const Polynomial & P1,
			const Polynomial & P2) const{

			res.resize(P1.size() + P2.size() - 1);

			for(int i = 0; i < res.size(); i++){
				_F.assign(res[i],_F.zero);
			}

			for(int i = 0; i < P1.size(); i++){
				for(int j = 0; j < P2.size(); j++){
					_F.axpyin(res[i + j],P1[i],P2[j]);
				}
			}

			return res;
		}
		//@}

		template<class Matrix1, class Matrix2>
		bool areEqual(const Matrix1 & A, const Matrix2 & B){
			if((A.rowdim() != B.rowdim()) || (A.coldim() != B.coldim())){
				return false ;
			}

			for(size_t i = 0 ; i < A.rowdim() ; ++i){
				for(size_t j = 0 ; j < A.coldim() ; ++j){
					if(!_F.areEqual(A.getEntry(i,j),B.getEntry(i,j))){ //!@bug use refs
						return false ;
					}
				}
			}

			return true ;
		}

		template<class Matrix>
		void setIdentity(Matrix & I){
			for(size_t i = 0 ; i < I.rowdim() ; ++i){
				for(size_t j = 0 ; j < I.coldim() ; ++j){
					if(i == j){
						I.setEntry(i,j,_F.one);
					}
					else{
						I.setEntry(i,j,_F.zero);
					}
				}
			}
		}

		template<class Matrix>
		void setZero(Matrix & I){
			// use Iterator
			for(size_t i = 0 ; i < I.rowdim() ; ++i){
				for(size_t j = 0 ; j < I.coldim() ; ++j){
					I.setEntry(i,j,_F.zero);
				}
			}
		}

		template<class Matrix1>
		bool isZero(const Matrix1 & A){
			for(size_t i = 0 ; i < A.rowdim() ; ++i){
				for(size_t j = 0 ; j < A.coldim() ; ++j){
					if(!_F.isZero(A.getEntry(i,j))){ //!@bug use refs
						return false;
					}
				}
			}

			return true ;
		}

		template<class Matrix1>
		bool isIdentity(const Matrix1 & A){
			if(A.rowdim() != A.coldim()){
				return false;
			}

			for(size_t i = 0 ; i < A.rowdim() ; ++i){
				if(!_F.isOne(A.getEntry(i,i))){
					return false;
				}
			}

			for(size_t i = 0 ; i < A.rowdim() ; ++i){
				for(size_t j = 0 ; j < i ; ++j){
					if(!_F.isZero(A.getEntry(i,j))){ //!@bug use refs
						return false;
					}
				}
			}

			for(size_t i = 0 ; i < A.rowdim() ; ++i){
				for(size_t j = i + 1 ; j < A.coldim() ; ++j){
					if(!_F.isZero(A.getEntry(i,j))){ //!@bug use refs
						return false;
					}
				}
			}

			return true ;
		}

		template<class Matrix1>
		bool isIdentityGeneralized(const Matrix1 & A){
			size_t mn = std::min(A.rowdim(),A.coldim());
			for(size_t i = 0 ; i < mn ; ++i){
				if(!_F.isOne(A.getEntry(i,i))){
					return false;
				}
			}

			for(size_t i = 0 ; i < A.rowdim() ; ++i){
				for(size_t j = 0 ; j < std::min(i,mn) ; ++j){
					if(!_F.isZero(A.getEntry(i,j))){ //!@bug use refs
						return false;
					}
				}
			}

			for(size_t i = 0 ; i < A.rowdim() ; ++i){
				for(size_t j = i+1 ; j < A.coldim() ; ++j){
					if(!_F.isZero(A.getEntry(i,j))){ //!@bug use refs
						return false;
					}
				}
			}

			return true;
		}

	//public:

		/** Print matrix.
		 * @param  os  Output stream to which matrix is written.
		 * @param  A   Matrix.
		 * @returns reference to os.
		 */
		template <class Matrix>
		inline std::ostream &write(std::ostream &os, const Matrix &A) const{
			return A.write(os, _F);
		}

		template <class Matrix>
		inline std::ostream &write(std::ostream &os,
		                           const Matrix &A,
		                           bool maple_format) const{

			return A.write(os, _F, maple_format);
		}

		/** Read matrix
		 * @param  is  Input stream from which matrix is read.
		 * @param  A   Matrix.
		 * @returns reference to is.
		 */
		template <class Matrix>
		inline std::istream &read(std::istream &is, Matrix &A) const{
			return A.read (is, _F);
		}

	}; /* end of class OpenCLMatrixDomain */

} /* end of namespace LinBox */

#ifdef __LINBOX_HAVE_OCL
	#include "linbox/algorithms/opencl-domain-factory.h"
	#include "linbox/algorithms/opencl-domain-util.inl"
	#include "linbox/algorithms/opencl-domain-memory.inl"
	#include "linbox/algorithms/opencl-domain.inl"
#endif

#endif // __LINBOX_opencl_matrix_domain_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
