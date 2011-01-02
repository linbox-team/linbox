/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/matrix/blas-matrix.h
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

/*! @file matrix/blas-matrix.h
 * @ingroup matrix
 * A \c BlasMatrix<\c _element > represents a matrix as an array of
 * <code>_element</code>s.
 *
 * @todo explain the differences with \c DenseMatrixBase<\c _element >
 * and \c DenseSubmatrix<\c _element > or
 * \c DenseMatrix<\c _field >...
 */

#ifndef __LINBOX_blas_matrix_H
#define __LINBOX_blas_matrix_H


#include <linbox/solutions/getentry.h>
#include <linbox/matrix/dense.h>
#include <linbox/matrix/dense-submatrix.h>
#include <linbox/util/debug.h>
#include <linbox/matrix/matrix-category.h>

#include "linbox/algorithms/cra-full-multip.h"
namespace LinBox
{

	// forward declaration
	template<class Element>
	class BlasMatrix;

	template <class Element>
	class MatrixContainerTrait<BlasMatrix<Element> > {
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};

	template <class Element>
	class MatrixContainerTrait<const BlasMatrix<Element> > {
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};

	/*! BlasMatrix.
	 * @ingroup matrix
	 * Limited docs so far.
	 */
	template <class _Element>
	class BlasMatrix : public DenseSubmatrix<_Element> {

	public:
		typedef _Element                                     Element;
		// typedef _Element* pointer ?
		typedef typename DenseSubmatrix<_Element>::Element * pointer;
		typedef BlasMatrix<Element>                           Self_t;

	protected:
		size_t   _stride;
		bool      _alloc;
		pointer     _ptr;

	private:

		/*! @internal
		 * @name create BlasMatrix
		 * @{ */

		/*! @internal
		 * Copy data according to blas container structure.
		 * Specialisation for BlasContainer.
		 */
		template <class _Matrix>
		void createBlasMatrix (const _Matrix& A, MatrixContainerCategory::BlasContainer)
		{
			typename _Matrix::ConstRawIterator         iter_value = A.rawBegin();
			RawIterator  iter_addr = this->rawBegin();
			for (;iter_value != A.rawEnd(); ++iter_value,++iter_addr)
				*iter_addr = *iter_value;
		}

		/*! @internal
		 * Copy data according to Matrix container structure.
		 * Specialisation for Container
		 */
		template <class Matrix>
		void createBlasMatrix (const Matrix& A, MatrixContainerCategory::Container)
		{
			// With both iterators, it is Segfaulting !!!!
			typename Matrix::ConstRawIndexedIterator  iter_index = A.rawIndexedBegin();
			for (;iter_index != A.rawIndexedEnd(); ++iter_index)
				this->_M->setEntry( iter_index.rowIndex(),
						    iter_index.colIndex(),
						    A.getEntry(iter_index.rowIndex(),iter_index.colIndex())
						  );
		}

		/*! @internal
		 * Copy data according to blackbox structure.
		 * Specialisation for Blackbox.
		 */
		template <class Matrix>
		void createBlasMatrix (const Matrix& A, MatrixContainerCategory::Blackbox)
		{
			typedef typename Matrix::Field Field;
			typename Field::Element one, zero;
			Field F = A.field();
			F. init(one, 1);
			F. init(zero, 0);

			std::vector<typename Field::Element> e(A.coldim(), zero), tmp(A.rowdim());
			typedef typename DenseSubmatrix<_Element>::ColIterator ColIterator ;
			ColIterator col_p;

			typename BlasMatrix<Element>::Col::iterator elt_p;

			typename std::vector<typename Field::Element>::iterator e_p, tmp_p;


			for (col_p = DenseSubmatrix<_Element>:: colBegin(), e_p = e.begin();
			     e_p != e.end(); ++ col_p, ++ e_p)
			{

				F.assign(*e_p, one);

				A.apply (tmp, e);

				for (tmp_p = tmp.begin(), elt_p = col_p -> begin();
				     tmp_p != tmp.end(); ++ tmp_p, ++ elt_p)

					F.assign(*elt_p, *tmp_p);

				F.assign(*e_p, zero);
			}
		}

		/*! @internal
		 * Copy data according to Matrix container structure (allow submatrix).
		 * Specialisation for Container
		 */
		template <class _Matrix>
		void createBlasMatrix (const _Matrix& A,
				       const size_t i0,const size_t j0,
				       const size_t m, const size_t n,
				       MatrixContainerCategory::Container)
		{

			typename _Matrix::ConstRawIterator         iter_value = A.rawBegin();
			typename _Matrix::ConstRawIndexedIterator  iter_index = A.rawIndexedBegin();

			for (;iter_value != A.rawEnd(); ++iter_value,++iter_index){
				size_t i,j;
				i=iter_index.rowIndex();
				j=iter_index.colIndex();
				if (( i >= i0) && (i< i0+m) && (j >= j0) && (j < j0+n))
					this->_M->setEntry(i-i0, j-j0, *iter_value);
			}
		}

		/*! @internal
		 * Copy data according to Matrix container structure (allow submatrix).
		 * Specialisation for BlasContainer.
		 */
		template <class Matrix>
		void createBlasMatrix (const Matrix& A,
				       const size_t i0,const size_t j0,
				       const size_t m, const size_t n,
				       MatrixContainerCategory::BlasContainer)
		{

			typename Matrix::ConstRawIterator         iter_value = A.rawBegin();
			typename Matrix::ConstRawIndexedIterator  iter_index = A.rawIndexedBegin();

			for (;iter_value != A.rawEnd(); ++iter_value,++iter_index){
				size_t i,j;
				i=iter_index.rowIndex();
				j=iter_index.colIndex();
				if (( i >= i0) && (i< i0+m) && (j >= j0) && (j < j0+n))
					this->_M->setEntry(i-i0, j-j0, *iter_value);
			}
		}

		/*! @internal
		 * Copy data according to blackbox structure (allow submatrix).
		 * Specialisation for Blackbox matrices
		 */
		template <class Matrix>
		void createBlasMatrix (const Matrix& A,
				       const size_t i0,const size_t j0,
				       const size_t m, const size_t n,
				       MatrixContainerCategory::Blackbox)
		{
			std::cerr << __func__ << ": not implemented yet" << std::flush << std::endl;
			exit(-1) ;
			//! @todo need to be implemented by succesive apply
			//! @todo lancer une exception générique "not implemented yet"
		}

		/*! @internal
		 * @}
		 */

	public:
		typedef typename DenseSubmatrix<_Element>::RawIterator RawIterator ; // for nvcc

		/* Constructors. */

		/*! Allocates a new \f$ 0 \times 0\f$ matrix.
		 */
		BlasMatrix () :
			DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (0,0)),0,0,0,0),
			_stride(0),  _alloc(true)
		{
			_ptr = this->_M->FullIterator();
		}

		/*! Allocates a new \f$ m \times n\f$ matrix.
		 * @param m rows
		 * @param n cols
		 */
		BlasMatrix (int m, int n) :
			DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (m,n)),0,0,m,n),
			_stride(n), _alloc(true)
		{
			_ptr = this->_M->FullIterator();
		}

		/*! Allocates a new \f$ m \times n\f$ matrix.
		 * @param m rows
		 * @param n cols
		 */
		BlasMatrix (size_t m, size_t n) :
			DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (m,n)),0,0,m,n),
			_stride(n), _alloc(true)
		{
			_ptr = this->_M->FullIterator();
		}

		/*! Constructor from a matrix stream.
		 * @param ms matrix stream.
		 */
		template< class Field >
		BlasMatrix(MatrixStream<Field>& ms) :
			DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (ms))),  _alloc(true)
		{
			_ptr = this->_M->FullIterator();
			_stride= this->coldim();
		}

		/*! Generic copy constructor from either a blackbox or a matrix container.
		 * @param A matrix to be copied
		 */
		template <class Matrix>
		BlasMatrix (const Matrix &A) :
			DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (A.rowdim(),A.coldim())),0,0,A.rowdim(),A.coldim()),
			_stride(A.coldim()) , _alloc(true)
		{
			_ptr = this->_M->FullIterator();
			createBlasMatrix(A, typename MatrixContainerTrait<Matrix>::Type());
		}

		/*! Generic copy constructor from either a blackbox or a matrix container (allow submatrix).
		 * @param A matrix to be copied
		 * @param i0
		 * @param j0
		 * @param m rows
		 * @param n columns
		 */
		template <class Matrix>
		BlasMatrix (const Matrix& A,
			    const size_t i0, const size_t j0,
			    const size_t m, const size_t n) :
			DenseSubmatrix<Element>( *(new DenseMatrixBase<Element> (A.rowdim(),A.coldim())),0,0,A.rowdim(),A.coldim()),
			_stride(A.coldim()) , _alloc(true)
		{
			_ptr = this->_M->FullIterator();
			createBlasMatrix(A, i0, j0, m, n, typename MatrixContainerTrait<Matrix>::Type());
		}

		/*! Constructor.
		 * @param A matrix to be copied
		 * @param F ?
		*/
		template<class _Matrix, class _Field>
		BlasMatrix (const _Matrix &A,  const _Field &F) :
			DenseSubmatrix<Element>( *(new DenseMatrixBase<Element> (A.rowdim(),A.coldim())),0,0,A.rowdim(),A.coldim() ),
			_stride(A.coldim()) ,
			_alloc(true)
		{
			_ptr = this->_M->FullIterator() ;
			rebind<_Field>()(*this,A,F);

		}

#if 1
		/*! Constructor from matrix (no copy).
		 * @param A DenseMatrixBase
		 */
		BlasMatrix ( DenseMatrixBase<Element>& A ) :
			DenseSubmatrix<Element>(A,0,0,A.rowdim(),A.coldim()),
			_stride(A.coldim()) , _alloc(false)
		{
			_ptr = this->_M->FullIterator();
		}

		/*! Constructor from matrix (no copy).
		 * @param A DenseMatrixBase
		 * @param i0
		 * @param j0
		 * @param m rows
		 * @param n columns
		 */
		BlasMatrix ( DenseMatrixBase<Element>& A,
			     const size_t i0, const size_t j0,
			     const size_t m, const size_t n) :
			DenseSubmatrix<Element>(A,i0,j0,m,n),
			_stride(A.coldim()) , _alloc(false)
		{
			_ptr = this->_M->FullIterator();
		}
#endif


		/*! Copy Constructor of a matrix (copying data).
		 * @param A matrix to be copied.
		 */
		BlasMatrix (const BlasMatrix<Element>& A) :
			DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (*A._M)),0,0,A.rowdim(),A.coldim()),
			_stride(A._stride), _alloc(true)
		{
			_ptr = this->_M->FullIterator();
		}


#if 1
		/// Copy Contructor of a matrix (no copy is done, just through pointer)
		/*! @param A BlasMatrix to be copied
		 */
		BlasMatrix(BlasMatrix<Element>& A) :
			DenseSubmatrix<Element>(A),
			_stride(A._stride), _alloc(false), _ptr(A._ptr)
		{}


		/// Copy Contructor of a submatrix (no copy is done, just through pointer)
		/*! @param A BlasMatrix to be copied
		 * @param i0
		 * @param j0
		 * @param m rows
		 * @param n columns
		 */
		BlasMatrix(BlasMatrix<Element>& A, const size_t i, const size_t j, const size_t m, const size_t n) :
			DenseSubmatrix<Element>(A,i,j,m,n),
			_stride(A._stride), _alloc(false),
			_ptr(A._ptr+ i*A._stride+j)
		{}
#endif



#if 0 /*  Create a blas matrix from a pointer of elements
       *  (without allocating a vector in DenseMatrixBase...)
       *  is it possible ? */
		BlasMatrix(size_t m, size_t n, size_t stride,
			   Element * A, size_t lda)
		{

		}
#endif


		/// Destructor.
		~BlasMatrix ()
		{
			if (_alloc)
				delete this->_M;
		}

		//! Rebind operator
		template<typename _Tp1>
		struct rebind {
			typedef BlasMatrix<typename _Tp1::Element> other;

			void operator() (other & Ap, const Self_t& A, const _Tp1& F)
			{

				typedef typename DenseSubmatrix<_Element>::ConstRawIndexedIterator ConstRawIndexedIterator ;
				typedef typename DenseSubmatrix<_Element>::ConstRawIterator ConstRawIterator ;
				ConstRawIterator         iter_value = A.rawBegin();
				ConstRawIndexedIterator  iter_index = A.rawIndexedBegin();
				typename _Tp1::Element tmp;
				for (;iter_value != A.rawEnd(); ++iter_value,++iter_index){
					F.init(  tmp, *iter_value );
					Ap.setEntry(iter_index.rowIndex(), iter_index.colIndex(),tmp);
				}
			}
		};


		//! operator = (copying data)
		BlasMatrix<Element>& operator= (const BlasMatrix<Element>& A)
		{

			DenseMatrixBase<Element> *tmp= this->_M;
			this->_M       = new DenseMatrixBase<Element>(*A._M);
			if (_alloc) {
				delete tmp;
			}
			this->_beg_row = A._beg_row;
			this->_end_row = A._end_row;
			this->_beg_col = A._beg_col;
			this->_end_col = A._end_col;
			_ptr     = this->_M->FullIterator();
			_alloc   = true;
			_stride  = A._stride;

			return *this;
		}

		/*! @internal
		 * Get read-only pointer to the matrix data.
		 */
		pointer getPointer() const
		{
			return _ptr;
		}

		/*! @internal
		 * Get write pointer to the matrix data.
		 * Data may be changed this way.
		 */
		pointer& getWritePointer()
		{
			return _ptr;
		}

		/*! @internal
		 * Get the stride of the matrix
		 */
		size_t getStride() const
		{
			return _stride;
		}

		/*! @internal
		 * Get a reference to the stride of the matrix.
		 * Modify stride this way.
		 */
		size_t& getWriteStride()
		{
			return _stride;
		}

		/*! @internal Is the matrix allocated ?
		 */
		bool isAllocated()
		{
			return _alloc ;
		}

	}; // end of class BlasMatrix


	/*! TAG for triangular blas matrix.
	* @see \ref FFLAS::FFLAS_DIAG and \ref FFLAS::FFLAS_UPLO enums in \c fflas/fflas.h.
	*/
	class BlasTag {
	public:
		typedef enum{low,up}       uplo; //!< upper or lower triangular
		typedef enum{unit,nonunit} diag; //!< unit or non unit diagonal
	};


	//! Triangular BLAS matrix.
	template <class Element>
	class TriangularBlasMatrix: public BlasMatrix<Element> {

	protected:

		BlasTag::uplo           _uplo; //!< upper or lower triangular
		BlasTag::diag           _diag; //!< unit or non unit diagonal

	public:

		/*! Constructor for a new \c TriangularBlasMatrix.
		 * @param m rows
		 * @param n cols
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		TriangularBlasMatrix (const size_t m, const size_t n,
				      BlasTag::uplo x=BlasTag::up, BlasTag::diag y= BlasTag::nonunit) :
			BlasMatrix<Element>(m, n ) , _uplo(x), _diag(y)
		{}

		/*! Constructor from a \c BlasMatrix (copy).
		 * @param A matrix
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		TriangularBlasMatrix (const BlasMatrix<Element>& A,
				      BlasTag::uplo x=BlasTag::up, BlasTag::diag y= BlasTag::nonunit) :
			BlasMatrix<Element>(A) , _uplo(x), _diag(y)
		{}

		/*! Constructor from a \c BlasMatrix (no copy).
		 * @param A matrix
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		TriangularBlasMatrix (BlasMatrix<Element>& A, BlasTag::uplo x=BlasTag::up, BlasTag::diag y= BlasTag::nonunit) :
			BlasMatrix<Element>(A), _uplo(x), _diag(y)
		{}

		/*! Constructor from a \c TriangularBlasMatrix (copy).
		 * @param A matrix
		 */
		TriangularBlasMatrix (const TriangularBlasMatrix<Element>& A) :
			BlasMatrix<Element>(A.rowdim(),A.coldim()), _uplo(A._uplo), _diag(A._diag)
		{
			switch (A._uplo) {
			case BlasTag::up:
				{
					for (size_t i=0;i<A.rowdim();++i)
						for (size_t j=i;j<A.coldim();++j)
							this->setEntry(i,j,A.getEntry(i,j));
					break;
				}
			case BlasTag::low:
				{
					for (size_t i=0;i<A.rowdim();++i) {
						for (size_t j=0;j<=i;++j)
							this->setEntry(i,j,A.getEntry(i,j));
					}

					break;
				}
			default:
				throw LinboxError ("Error in copy constructor of TriangularBlasMatrix (incorrect argument)");
			}
		}

		/*! Generic constructor from a \c Matrix (no copy).
		 * @param A matrix
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		template<class Matrix>
		TriangularBlasMatrix (const Matrix& A, BlasTag::uplo x=BlasTag::up, BlasTag::diag y= BlasTag::nonunit) :
			BlasMatrix<Element>(A.rowdim(),A.coldim()), _uplo(x), _diag(y)
		{
			switch (x) {
			case BlasTag::up:
				{
					for (size_t i=0;i<A.rowdim();++i){
						for (size_t j=i;j<A.coldim();++j) {
							Element tmp;
							this->setEntry(i,j,getEntry(tmp, A,i,j));
						}
					}
					break;
				}
			case BlasTag::low:
				{
					for (size_t i=0;i<A.rowdim();++i) {
						for (size_t j=0;j<=i;++j) {
							Element tmp;
							this->setEntry(i,j,getEntry(tmp, A,i,j));
						}
					}

					break;
				}
			default:
				throw LinboxError ("Error in copy constructor of TriangularBlasMatrix (incorrect argument)");
			}
		}

		/// get the shape of the matrix (upper or lower)
		BlasTag::uplo getUpLo() const { return _uplo;}

		/// Is the diagonal implicitly unit ?
		BlasTag::diag getDiag() const { return _diag;}

	}; // end of class TriangularBlasMatrix

	template <class Element>
	struct MatrixTraits< BlasMatrix<Element> > {
		typedef BlasMatrix<Element> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};

	template <class Element>
	struct MatrixTraits< const BlasMatrix<Element> > {
		typedef const BlasMatrix<Element> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};


	/** Class used for permuting indices.
	 * For example, create a vector <code>(0 1 2 ...)</code> over \c
	 * size_t, then apply a permutation to it using a \c BlasMatrixDomain to
	 * get the natural representation of the permutation.
	 */
	class indexDomain {
	public:
		typedef size_t Element;
	public:
		indexDomain() {};
		template <class ANY>
		size_t init(size_t& dst, const ANY& src) const {
			return dst = static_cast<size_t>(src);
		}
		template <class ANY>
		size_t assign(ANY& dst, const size_t& src) const {
			return dst = static_cast<ANY>(src);
		}
	};

#if 0 /*  original BlasPermutation */
	// Dan Roche 7-8-04 Changed _P to _PP to avoid confict with a macro defined in
	// <iostream> somewhere.
	/*! BlasPermutation.
	 * Lapack ipiv style compressed permutation.
	 * @todo to be deprecated
	 */
	class BlasPermutation {


	public:

		/// null constructor
		BlasPermutation() {};

		BlasPermutation( const size_t n ) :
			_PP(n), _order( n ), _size( -1 )
		{};

		BlasPermutation( const std::vector<size_t> P ) :
			_PP( P ), _order( P.size() ) , _size( -1 )
		{};

		BlasPermutation( const std::vector<size_t> P, size_t order ) :
			_PP( P ), _order( order ) , _size( -1 )
		{};

		BlasPermutation( const std::vector<size_t> P,
				 size_t order, size_t size ) :
			_PP( P ), _order( order ) , _size( size )
		{};

		BlasPermutation( const BlasPermutation& P) :
			_PP( P._PP ), _order( P._order ) , _size( size )
		{};

		BlasPermutation& operator=( const BlasPermutation& P )
		{
			_PP    = P._PP;
			_order = P._order;
			_size  = P._size ;
			return *this;
		}

		const size_t* getPointer() const
		{
			return &_PP[0];
		}

		size_t* getWritePointer()
		{
			return &_PP[0];
		}

		size_t  getOrder()  const
		{
			return _order;
		}

		size_t  getSize()  const
		{
			if (_size == -1)
				_size = (*(std::max_element(_PP.begin(),_PP.end())))+1 ;
			return _size;
		}

		BlasPermutation& extendTrivially(const size_t newSize)
		{
			if (newSize < _order)
				std::cerr << "WARNING: attempting to reduce size of permutation.";
			_PP.resize(newSize);
			for (size_t i=_order; i<newSize; i++)
				_PP[i] = i;
			_order = newSize;
			return *this;
		};

	protected:

		std::vector<size_t>  _PP ;    //!< Lapack representation of the permutation
		size_t               _order;  //!< size of the representation (number of permutations)
		mutable size_t       _size ;  //!< size of the permutation (computed if necessary)

	private :

		/// compresses PackedPermutation to a smaller \c r_.
		void Compress_()
		{
			if (_order==0) {
				_PP.resize(0) ;
				// Id_ = true ;
				return ;
			}
			linbox_check(_order);
			size_t rr = _order-1 ;
			while ( rr && (_PP[rr] == 0  )) --rr ;    // removing trailing zeros
			while ( rr && (_PP[rr] == rr )) --rr ;    // useless information
			if ((rr == 0) && (_PP[0] == 0)) {
				_order = 0 ;
				_size = 0  ;
				Id_ = true ;
				_PP.resize(0) ;
				return ;
			}
			_order = rr+1 ;
			_PP.resize(_order,0);   // done cleaning.
			// recomputing n_ if lost.
			if (_size !=  -1) {
				_size = getSize();
			}
			cleaned_ = true ;
			return ;
		}


	}; // end of class BlasPermutation
#endif

	/*! TransposedBlasMatrix.
	 * NO DOC
	 */
	template< class Matrix >
	class TransposedBlasMatrix {

	public:

		/*! NO DOC
		 * @param M
		 */
		TransposedBlasMatrix ( Matrix& M ) :
			_M(M)
		{}

		/*! NO DOC
		 */
		Matrix& getMatrix() const { return _M; }

	protected:
		Matrix& _M; //!< NO DOC
	};

	/*! TransposedBlasMatrix.
	 * NO DOC
	 */
#if !defined(__INTEL_COMPILER) && !defined(__CUDACC__)
	template <>
#endif
	template< class Matrix >
	class TransposedBlasMatrix< TransposedBlasMatrix< Matrix > > : public Matrix {

	public:
		/*! TransposedBlasMatrix.
		 * NO DOC
		 */
		TransposedBlasMatrix ( Matrix& M ) :
			Matrix(M)
		{}

		/*! TransposedBlasMatrix.
		 * NO DOC
		 */
		TransposedBlasMatrix ( const Matrix& M ) :
			Matrix(M)
		{}

	};


} // end of namespace LinBox

#endif // __LINBOX_blas_matrix_H

