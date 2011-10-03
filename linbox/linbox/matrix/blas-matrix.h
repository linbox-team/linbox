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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/*! @file matrix/blas-matrix.h
 * @ingroup matrix
 * A \c BlasMatrix<\c _element > represents a matrix as an array of
 * <code>_element</code>s.
 *
 */

#ifndef __LINBOX_blas_matrix_H
#define __LINBOX_blas_matrix_H

#include "linbox/vector/subiterator.h"
#include "linbox/vector/subvector.h"

#include "linbox/algorithms/linbox-tags.h"
#include "linbox/util/debug.h"
#include "linbox/matrix/matrix-category.h"

// Blas Matrix
namespace LinBox
{
	template<class _Element>
	class BlasSubmatrix ;

	/*! Dense matrix representation.
	 * @ingroup matrix
	 * A BlasMatrix is a matrix of \p _Element, with the structure of BLAS matrices.
	 * It is basically a vector of \p _Element.
	 * In the Mother model, a BlasMatrix is allocated by the user.
	 */
	template <class _Element>
	class BlasMatrix {
		// private :

	public:
		typedef _Element                            Element;    //!< Element type
		typedef typename RawVector<Element>::Dense      Rep;    //!< Actually a <code>std::vector<Element></code> (or alike.)
		// typedef _Element *                          pointer;    //!< pointer to elements type
		typedef typename Rep::pointer               pointer;    //!< pointer type to elements
		typedef const pointer                 const_pointer;    //!< const pointer type
		typedef BlasMatrix<Element>                  Self_t;    //!< Self type

	protected:
		size_t      _row;
		size_t      _col;
		Rep         _rep;
		pointer     _ptr;

	private:

#if 0
		void makePointer()
		{
#if 0
			if (_row && _col) {
				_ptr = malloc( _row*_col*sizeof(_Element) ) ;
				linbox_check(_ptr);
			}
			else
				_ptr = NULL ;
#endif
			_rep = Rep(_row*_col);
			_ptr = &_rep[0];
		}
#endif

		/*! @internal
		 * @name create BlasMatrix
		 * @{ */

		/*! @internal
		 * Copy data according to blas container structure.
		 * Specialisation for BlasContainer.
		 */
		void createBlasMatrix (const BlasMatrix<Element> & A) ;

		/*! @internal
		 * Copy data according to blas container structure.
		 * Specialisation for BlasContainer.
		 */
		template <class _Matrix>
		void createBlasMatrix (const _Matrix& A, MatrixContainerCategory::BlasContainer) ;

		/*! @internal
		 * Copy data according to Matrix container structure.
		 * Specialisation for Container
		 */
		template <class Matrix>
		void createBlasMatrix (const Matrix& A, MatrixContainerCategory::Container) ;

		/*! @internal
		 * Copy data according to blackbox structure.
		 * Specialisation for Blackbox.
		 */
		template <class Matrix>
		void createBlasMatrix (const Matrix& A, MatrixContainerCategory::Blackbox) ;

		/*! @internal
		 * Copy data according to Matrix container structure (allow submatrix).
		 * Specialisation for Container
		 */
		template <class _Matrix>
		void createBlasMatrix (const _Matrix& A,
				       const size_t i0,const size_t j0,
				       const size_t m, const size_t n,
				       MatrixContainerCategory::Container) ;

		/*! @internal
		 * Copy data according to Matrix container structure (allow submatrix).
		 * Specialisation for BlasContainer.
		 */
		template <class Matrix>
		void createBlasMatrix (const Matrix& A,
				       const size_t i0,const size_t j0,
				       const size_t m, const size_t n,
				       MatrixContainerCategory::BlasContainer) ;

		/*! @internal
		 * Copy data according to blackbox structure (allow submatrix).
		 * Specialisation for Blackbox matrices
		 * @todo need to be implemented by succesive apply
		 */
		template <class Matrix>
		void createBlasMatrix (const Matrix& A,
				       const size_t i0,const size_t j0,
				       const size_t m, const size_t n,
				       MatrixContainerCategory::Blackbox) ;

		void createBlasMatrix (const Element * v) ;
		void createBlasMatrix (const std::vector<Element> & v) ;
		/*! @internal
		 * @}
		 */

	public:

		//////////////////
		// CONSTRUCTORS //
		//////////////////


		/*! Allocates a new \f$ 0 \times 0\f$ matrix.
		*/
		BlasMatrix () ;

		/*! Allocates a new \f$ m \times n\f$ matrix.
		 * @param m rows
		 * @param n cols
		 */
		//@{

		template<class T>
		BlasMatrix ( unsigned int m, T n) ;

		template<class T>
		BlasMatrix (long m, T n) ;

		template<class T>
		BlasMatrix (unsigned long m, T  n) ;

		template<class T>
		BlasMatrix (int m, T n) ;

		template<class T>
		BlasMatrix ( Integer & m, T n) ;

		//@}


		/*! Constructor from a matrix stream.
		 * @param ms matrix stream.
		 */
		template< class Field >
		BlasMatrix(MatrixStream<Field>& ms) ;

		/*! Generic copy constructor from either a blackbox or a matrix container.
		 * @param A matrix to be copied
		 */
		template <class Matrix>
		BlasMatrix (const Matrix &A) ;

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
			    const size_t m,  const size_t n) ;

		/*! Constructor.
		 * @param A matrix to be copied
		 * @param F Field
		 */
		template<class _Matrix, class _Field>
		BlasMatrix (const _Matrix &A,  const _Field &F) ;

		/*! Copy Constructor of a matrix (copying data).
		 * @param A matrix to be copied.
		 */
		BlasMatrix (const BlasMatrix<Element>& A) ;

		/*- Copy Constructor of a matrix (copying data).
		 * @param A matrix to be copied.
		 */
		// BlasMatrix (const BlasSubmatrix<Element>& A) ;

		/*! Create a BlasMatrix from a vector of elements
		 * @param A matrix to be copied.
		 */
		BlasMatrix (const std::vector<Element>& v, size_t m, size_t n) ;

		/*! Create a BlasMatrix from an array of elements
		 * @param A matrix to be copied.
		 */
		BlasMatrix (const _Element * v, size_t m, size_t n) ;

		/// Destructor.
		~BlasMatrix () ;

		//! operator = (copying data)
		BlasMatrix<Element>& operator= (const BlasMatrix<Element>& A) ;

		//! Rebind operator
		template<typename _Tp1>
		struct rebind ;

		//////////////////
		//  DIMENSIONS  //
		//////////////////

		/** Get the number of rows in the matrix.
		 * @returns Number of rows in matrix
		 */
		size_t rowdim() const ;

		/** Get the number of columns in the matrix.
		 * @returns Number of columns in matrix
		 */
		size_t coldim() const ;

		/*! Get the stride of the matrix.
		 */
		size_t getStride() const;

		/*!Get a reference to the stride of the matrix.
		 * Modify stride this way.
		 */
		size_t& getWriteStride();


		/** Resize the matrix to the given dimensions.
		 * The state of the matrix's entries after a call to this method is
		 * undefined
		 * @param m Number of rows
		 * @param n Number of columns
		 * @param val
		 */
		void resize (size_t m, size_t n, const Element& val = Element()) ;

		//////////////////
		//   ELEMENTS   //
		//////////////////

		/*! @internal
		 * Get read-only pointer to the matrix data.
		 */
		pointer getPointer() const ;

		const_pointer &getConstPointer() const ;


		/*! @internal
		 * Get write pointer to the matrix data.
		 * Data may be changed this way.
		 */
		pointer& getWritePointer() ;

		/** Set the entry at the (i, j) position to a_ij.
		 * @param i Row number, 0...rowdim () - 1
		 * @param j Column number 0...coldim () - 1
		 * @param a_ij Element to set
		 */
		void setEntry (size_t i, size_t j, const Element &a_ij) ;

		/** Get a writeable reference to the entry in the (i, j) position.
		 * @param i Row index of entry
		 * @param j Column index of entry
		 * @returns Reference to matrix entry
		 */
		Element &refEntry (size_t i, size_t j) ;

		/** Get a read-only reference to the entry in the (i, j) position.
		 * @param i Row index
		 * @param j Column index
		 * @returns Const reference to matrix entry
		 */
		const Element &getEntry (size_t i, size_t j) const ;

		/** Copy the (i, j) entry into x, and return a reference to x.
		 * This form is more in the Linbox style and is provided for interface
		 * compatibility with other parts of the library
		 * @param x Element in which to store result
		 * @param i Row index
		 * @param j Column index
		 * @returns Reference to x
		 */
		Element &getEntry (Element &x, size_t i, size_t j) const ;

		///////////////////
		// TRANSPOSE &AL //
		///////////////////

		/*! Creates a transposed matrix of \c *this.
		 * @param[in] tM
		 * @return the transposed matrix of this.
		 */
		BlasMatrix<Element> transpose(BlasMatrix<Element> & tM) const ;


		/*! Transpose (inplace).
		 * If rows and columns agree, we can transpose inplace.
		 */
		void transpose() ;

		/*! Reverse the rows of a matrix.
		 * This is done inplace.
		 * Let J=antiDiag(1) (or the matrix of the reverse
		 * permutation or the matrix (i,j) = (i+j+1==m)). Then,
		 * we compute A <- J.A;
		 */
		void reverseRows() ;

		/*! Reverse the columns of a matrix.
		 * This is done inplace.
		 * This is A <- J.A
		 */
		void reverseCols() ;

		/*! Reverse the rows/columns of a matrix.
		 * This is done inplace.
		 * This is A <- J.A.J
		 */
		void reverse() ;

		///////////////////
		//      I/O      //
		///////////////////

		/** Read the matrix from an input stream.
		 * The stream is in SMS or DENSE format
		 * @param file Input stream from which to read
		 * @param F Field over which to read
		 */
		template <class Field>
		std::istream &read (std::istream &file, const Field &F);

		/** Write the matrix to an output stream.
		 * @todo factorise writing matrices code.
		 * @param os Output stream to which to write
		 * @param F Field over which to write
		 * @param mapleFormat write in Maple format ?
		 */
		template <class Field>
		std::ostream &write (std::ostream &os, const Field &F,
				     bool mapleFormat=true) const;

		/** Write brutally the matrix to an output stream.
		 * This a raw version of \c write(os,F) (no field is given).
		 * @param os Output stream to which to write
		 * @param mapleFormat write in maple format ?
		 */
		std::ostream &write (std::ostream &os,
				     bool mapleFormat=true) const;


		///////////////////
		//   ITERATORS   //
		///////////////////

		/** @name Column of rows iterator
		 * \brief
		 * The column of rows iterator traverses the rows of the
		 * matrix in ascending order. Dereferencing the iterator yields
		 * a row vector in dense format
		 */
		//@{
		typedef Subvector<typename Rep::iterator, typename Rep::const_iterator> Row;
		typedef Subvector<typename Rep::const_iterator>                    ConstRow;

		/*!  Row Iterator.
		 * @ingroup iterators
		 * @brief NO DOC
		 */
		class RowIterator;
		/*! Const Row Iterator.
		 * @ingroup iterators
		 * @brief NO DOC
		 */
		class ConstRowIterator;

		RowIterator      rowBegin ();
		RowIterator      rowEnd ();
		ConstRowIterator rowBegin () const;
		ConstRowIterator rowEnd   () const;
		//@}

		/** @name Row of columns iterator
		 * \brief
		 * The row of columns iterator traverses the columns of the
		 * matrix in ascending order. Dereferencing the iterator yields
		 * a column vector in dense format
		 */
		//@{
		typedef Subvector<Subiterator<typename Rep::iterator> >            Col;
		typedef Subvector<Subiterator<typename Rep::const_iterator> > ConstCol;
		typedef Col           Column;
		typedef ConstCol ConstColumn;

		/*! Col Iterator.
		 * @ingroup iterators
		 * @brief NO DOC
		 */
		class ColIterator;
		/*! Const Col Iterator.
		 * @ingroup iterators
		 * @brief NO DOC
		 */
		class ConstColIterator;

		ColIterator      colBegin ();
		ColIterator      colEnd ();
		ConstColIterator colBegin () const;
		ConstColIterator colEnd ()   const;
		//@}

		/** @name Iterator
		 * \brief
		 *
		 * The iterator is a method for accessing all entries in the matrix
		 * in some unspecified order. This can be used, e.g. to reduce all
		 * matrix entries modulo a prime before passing the matrix into an
		 * algorithm.
		 */
		//@{
		typedef typename Rep::iterator Iterator;
		typedef typename Rep::const_iterator ConstIterator;

		Iterator      Begin ();
		Iterator      End   ();
		ConstIterator Begin () const;
		ConstIterator End   () const;
		//@}

		/** @name Raw Indexed iterator
		 * \brief
		 *
		 * Like the raw iterator, the indexed iterator is a method for
		 * accessing all entries in the matrix in some unspecified order.
		 * At each position of the the indexed iterator, it also provides
		 * the row and column indices of the currently referenced entry.
		 * This is provided through it's \c rowIndex() and \c colIndex() functions.
		 */
		//@{
		class IndexedIterator;
		/*! Const Indexed Iterator.
		 * @ingroup iterators
		 * @brief NO DOC
		 */
		class ConstIndexedIterator;

		IndexedIterator      IndexedBegin ();
		IndexedIterator      IndexedEnd   ();
		ConstIndexedIterator IndexedBegin () const;
		ConstIndexedIterator IndexedEnd   () const;
		//@}

		/** Retrieve a reference to a row.
		 * Since rows may also be indexed, this allows A[i][j] notation
		 * to be used.
		 * @param i Row index
		 */
		//@{
		Row      operator[] (size_t i) ;
		ConstRow operator[] (size_t i) const ;
		//@}

		///////////////////
		//     MISC     //
		///////////////////


		/** Compute column density.
		 * @param v
		 */
		template <class Vector>
		Vector &columnDensity (Vector &v) const ;

	}; // end of class BlasMatrix

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

	/*! Write a matrix to a stream.
	 * The C++ way using <code>operator<<</code>
	 * @param o output stream
	 * @param M matrix to write.
	 */
	template<class T>
	std::ostream& operator<< (std::ostream & o, const BlasMatrix<T> & Mat)
	{
		return Mat.write(o);
	}

} // end of namespace LinBox

// Blas Submatrix
namespace LinBox
{
	/*! Dense Submatrix representation.
	 * @ingroup matrix
	 * A @ref BlasSubmatrix is a matrix of \p _Element, with the structure of BLAS matrices.
	 * It is basically a read/write view on a vector of \p _Element.
	 * In the Mother model, a @ref BlasSubmatrix is not allocated.
	 * <p>
	 * This matrix type conforms to the same interface as @ref BlasMatrix,
	 * except that you cannot resize it. It represents a submatrix of a dense
	 * matrix. Upon construction, one can freely manipulate the entries in the
	 * DenseSubmatrix, and the corresponding entries in the underlying
	 * @ref BlasMatrix will be modified.


	 */
	template <class _Element>
	class BlasSubmatrix {
	public :
		typedef _Element                  Element;      //!< Element type
		typedef BlasSubmatrix<_Element>   Self_t;       //!< Self type


	protected:
		BlasMatrix<Element> *_M;       //!< Parent BlasMatrix (ie raw vector)
		size_t _row;                   //!< row dimension of Submatrix
		size_t _col;                   //!< col dimension of Submatrix
		size_t _r0;                    //!< upper left corner row of Submatrix in \p _M
		size_t _c0;                    //!< upper left corner row of Submatrix in \p _M
		size_t _stride ;               //!< number of columns in \p _M (or stride of \p _M)
		size_t _off;

	public:

		//////////////////
		// CONSTRUCTORS //
		//////////////////


		/*  constructors */

		/** NULL constructor.  */
		BlasSubmatrix () ;

		/** Constructor from an existing @ref BlasMatrix and dimensions.
		 * \param M Pointer to @ref BlasMatrix of which to construct submatrix
		 * \param rowbeg Starting row
		 * \param colgeb Starting column
		 * \param Rowdim Row dimension
		 * \param Coldim Column dimension
		 */
		BlasSubmatrix (const BlasMatrix<Element> &M,
			       size_t rowbeg,
				size_t colbeg,
				size_t Rowdim,
				size_t Coldim);

		/** Constructor from an existing @ref DenseMatrixBase
		 * \param M Pointer to @ref DenseMatrixBase of which to construct submatrix
		 */
		BlasSubmatrix (const BlasMatrix<Element> &M);


		/** Constructor from an existing submatrix and dimensions
		 * @param SM Constant reference to BlasSubmatrix from which to
		 *           construct submatrix
		 * @param rowbeg Starting row
		 * @param colbeg Starting column
		 * @param Rowdim Row dimension
		 * @param Coldim Column dimension
		 */
		BlasSubmatrix (const BlasSubmatrix<Element> &SM,
				size_t rowbeg,
				size_t colbeg,
				size_t Rowdim,
				size_t Coldim);

		/** Copy constructor.
		 * @param SM Submatrix to copy
		 */
		BlasSubmatrix (const BlasSubmatrix<Element> &SM);


		/*  Members  */

		/** Assignment operator.
		 * Assign the given submatrix to this one
		 * This is <i>only</i> renaming !
		 * There is no copy because BlasSubmatrix owns nothing.
		 * @param SM Submatrix to assign
		 * @return Reference to this submatrix
		 */
		BlasSubmatrix &operator = (const BlasSubmatrix<Element> &SM);

		template<typename _Tp1>
		struct rebind {
			typedef BlasSubmatrix<typename _Tp1::Element> other;
		};

		//////////////////
		//  DIMENSIONS  //
		//////////////////

		/** Get the number of rows in the matrix
		 * @return Number of rows in matrix
		 */
		size_t rowdim () const;

		/** Get the number of columns in the matrix
		 * @return Number of columns in matrix
		 */
		size_t coldim () const ;

		/*! Get the stride of the matrix.
		 * @return stride of submatrix (number of cols of dense base matrix)
		 */
		size_t getStride() const;


		///////////////////
		//      I/O      //
		///////////////////

		/** Read the matrix from an input stream.
		 * @param file Input stream from which to read
		 * @param field
		 */
		template<class Field>
		std::istream& read (std::istream &file, const Field& field);

		/** Write the matrix to an output stream.
		 * @param os Output stream to which to write
		 * @param field
		 * @param mapleFormat write in Maple(r) format ?
		 */
		template<class Field>
		std::ostream& write (std::ostream &os, const Field& field,
				     bool mapleFormat = true) const;

		/** Write the matrix to an output stream.
		 * This a raw version of \c write(os,F) (no field is given).
		 * @param os Output stream to which to write
		 * @param mapleFormat write in Maple(r) format ?
		 */
		std::ostream& write (std::ostream &os,
				     bool mapleFormat = true) const;


		//////////////////
		//   ELEMENTS   //
		//////////////////


		/** Set the entry at (i, j).
		 * @param i Row number, 0...rowdim () - 1
		 * @param j Column number 0...coldim () - 1
		 * @param a_ij Element to set
		 */
		void setEntry (size_t i, size_t j, const Element &a_ij) ;

		/** Get a writeable reference to an entry in the matrix.
		 * @param i Row index of entry
		 * @param j Column index of entry
		 * @return Reference to matrix entry
		 */
		Element &refEntry (size_t i, size_t j) ;

		/** Get a read-only individual entry from the matrix.
		 * @param i Row index
		 * @param j Column index
		 * @return Const reference to matrix entry
		 */
		const Element &getEntry (size_t i, size_t j) const ;

		/** Get an entry and store it in the given value.
		 * This form is more in the Linbox style and is provided for interface
		 * compatibility with other parts of the library
		 * @param x Element in which to store result
		 * @param i Row index
		 * @param j Column index
		 * @return Reference to x
		 */
		Element &getEntry (Element &x, size_t i, size_t j) const ;


		///////////////////
		//   ITERATORS   //
		///////////////////

		//! @name Forward declaration of Raw Iterators.
		//@{
		class Iterator  ;
		class ConstIterator ;

		class IndexedIterator ;
		class ConstIndexedIterator ;
		//@}


		/** @name typedef'd Row Iterators.
		 *\brief
		 * The row iterator gives the rows of the
		 * matrix in ascending order. Dereferencing the iterator yields
		 * a row vector in dense format
		 * @{
		 */
		typedef typename BlasMatrix<Element>::RowIterator            RowIterator;
		typedef typename BlasMatrix<Element>::ConstRowIterator       ConstRowIterator;
		typedef typename BlasMatrix<Element>::Row                    Row;
		typedef typename BlasMatrix<Element>::ConstRow               ConstRow;
		//@} Row Iterators

		/** @name typedef'd Column Iterators.
		 *\brief
		 * The columns iterator gives the columns of the
		 * matrix in ascending order. Dereferencing the iterator yields
		 * a column vector in dense format
		 * @{
		 */
		typedef typename BlasMatrix<Element>::ColIterator            ColIterator;
		typedef typename BlasMatrix<Element>::ConstColIterator       ConstColIterator;
		typedef typename BlasMatrix<Element>::Col                    Col;
		typedef typename BlasMatrix<Element>::Column                 Column;
		typedef typename BlasMatrix<Element>::ConstCol               ConstCol;
		//@} // Column Iterators



		RowIterator      rowBegin ();        //!< iterator to the begining of a row
		RowIterator      rowEnd ();          //!< iterator to the end of a row
		ConstRowIterator rowBegin () const;  //!< const iterator to the begining of a row
		ConstRowIterator rowEnd ()   const;  //!< const iterator to the end of a row

		ColIterator      colBegin ();
		ColIterator      colEnd ();
		ConstColIterator colBegin () const;
		ConstColIterator colEnd ()   const;

		Iterator      Begin ();
		Iterator      End ();
		ConstIterator Begin () const;
		ConstIterator End ()   const;


		IndexedIterator      IndexedBegin();
		IndexedIterator      IndexedEnd();
		ConstIndexedIterator IndexedBegin() const;
		ConstIndexedIterator IndexedEnd()   const;

		/*!  operator[].
		 * Retrieve a reference to a row
		 * @param i Row index
		 */
		Row      operator[] (size_t i) ;
		ConstRow operator[] (size_t i) const ;

	};


	template <class Element>
	struct MatrixTraits< BlasSubmatrix<Element> > {
		typedef BlasSubmatrix<Element> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};

	template <class Element>
	struct MatrixTraits< const BlasSubmatrix<Element> > {
		typedef const BlasSubmatrix<Element> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};

	/*! Write a matrix to a stream.
	 * The C++ way using <code>operator<<</code>
	 * @param o output stream
	 * @param M matrix to write.
	 */
	template<class T>
	std::ostream& operator<< (std::ostream & o, const BlasSubmatrix<T> & Mat)
	{
		return Mat.write(o);
	}

}

// Triangular, Transposed Matrix
namespace LinBox
{
	//! Triangular BLAS matrix.
	template <class Element>
	class TriangularBlasMatrix: public BlasMatrix<Element> {

	protected:

		LinBoxTag::Shape          _uplo; //!< upper or lower triangular
		LinBoxTag::Diag           _diag; //!< unit or non unit diagonal

	public:

		/*! Constructor for a new \c TriangularBlasMatrix.
		 * @param m rows
		 * @param n cols
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		TriangularBlasMatrix (const size_t m, const size_t n,
				      LinBoxTag::Shape x=LinBoxTag::Upper, LinBoxTag::Diag y= LinBoxTag::NonUnit) :
			BlasMatrix<Element>(m, n ) , _uplo(x), _diag(y)
		{}

		/*! Constructor from a \c BlasMatrix (copy).
		 * @param A matrix
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		TriangularBlasMatrix (const BlasMatrix<Element>& A,
				      LinBoxTag::Shape x=LinBoxTag::Upper, LinBoxTag::Diag y= LinBoxTag::NonUnit) :
			BlasMatrix<Element>(A) , _uplo(x), _diag(y)
		{}

		/*! Constructor from a \c BlasMatrix (no copy).
		 * @param A matrix
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		TriangularBlasMatrix (BlasMatrix<Element>& A, LinBoxTag::Shape x=LinBoxTag::Upper, LinBoxTag::Diag y= LinBoxTag::NonUnit) :
			BlasMatrix<Element>(A), _uplo(x), _diag(y)
		{}

		/*! Constructor from a \c TriangularBlasMatrix (copy).
		 * @param A matrix
		 */
		TriangularBlasMatrix (const TriangularBlasMatrix<Element>& A) :
			BlasMatrix<Element>(A.rowdim(),A.coldim()), _uplo(A._uplo), _diag(A._diag)
		{
			switch (A._uplo) {
			case LinBoxTag::Upper:
				{
					for (size_t i=0;i<A.rowdim();++i)
						for (size_t j=i;j<A.coldim();++j)
							this->setEntry(i,j,A.getEntry(i,j));
					break;
				}
			case LinBoxTag::Lower:
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
		TriangularBlasMatrix (const Matrix& A, LinBoxTag::Shape x=LinBoxTag::Upper, LinBoxTag::Diag y= LinBoxTag::NonUnit) :
			BlasMatrix<Element>(A.rowdim(),A.coldim()), _uplo(x), _diag(y)
		{
			switch (x) {
			case LinBoxTag::Upper:
				{
					for (size_t i=0;i<A.rowdim();++i){
						for (size_t j=i;j<A.coldim();++j) {
							Element tmp;
							this->setEntry(i,j,getEntry(tmp, A,i,j));
						}
					}
					break;
				}
			case LinBoxTag::Lower:
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
		LinBoxTag::Shape getUpLo() const { return _uplo;}

		/// Is the diagonal implicitly unit ?
		LinBoxTag::Diag getDiag() const { return _diag;}

	}; // end of class TriangularBlasMatrix

} // LinBox

namespace LinBox
{


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
		int characteristic() const { return 0 ; }
	};
}


namespace LinBox
{
	/*! TransposedBlasMatrix.
	 * NO DOC
	 */
	template< class Matrix >
	class TransposedBlasMatrix {

	public:

		/*! NO DOC
		 * @param M
		 */
		TransposedBlasMatrix ( Matrix& Mat ) :
			_M(Mat)
		{}

		/*! NO DOC
		*/
		Matrix& getMatrix() const
		{
			return _M;
		}

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
		TransposedBlasMatrix ( Matrix& Mat ) :
			Matrix(Mat)
		{}

		/*! TransposedBlasMatrix.
		 * NO DOC
		 */
		TransposedBlasMatrix ( const Matrix& Mat ) :
			Matrix(Mat)
		{}

	};


}

#include "blas-matrix.inl"
#include "blas-submatrix.inl"

#endif // __LINBOX_blas_matrix_H

