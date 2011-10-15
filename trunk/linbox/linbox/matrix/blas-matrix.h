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
 * A \c BlasMatrix<\c _Field > represents a matrix as an array of
 * <code>_Field::Element</code>s. It also has the BlasBlackbox interface.
 *
 */

#ifndef __LINBOX_blas_matrix_H
#define __LINBOX_blas_matrix_H

#include "linbox/vector/subiterator.h"
#include "linbox/vector/subvector.h"

#include "linbox/util/debug.h"
#include "linbox/matrix/matrix-category.h"
#include "linbox/algorithms/linbox-tags.h"

namespace LinBox
{
	namespace Protected
	{

		template <class Field>
		bool checkBlasApply(const Field &F, size_t n)
		{

			integer chara, card;
			F.characteristic(chara);
			F.cardinality(card);

			if ((chara != card) || chara == 0)
				return false;
			else
				if (n*chara*chara < integer("9007199254740992"))
					return true;
				else
					return false;
		}
	}
}


// Blas Matrix
namespace LinBox
{
	template<class _Field>
	class BlasSubmatrix ;

	/*! Dense matrix representation.
	 * @ingroup matrix
	 * A BlasMatrix is a matrix of \p _Field::Element, with the structure of BLAS matrices.
	 * It is basically a vector of \p _Field::Element.
	 * In the Mother model, a BlasMatrix is allocated by the user.
	 */
	template <class _Field>
	class BlasMatrix {
		// private :

	public:
		typedef _Field                                Field;
		typedef typename Field::Element             Element;    //!< Element type
		typedef typename RawVector<Element>::Dense      Rep;    //!< Actually a <code>std::vector<Element></code> (or alike.)
		typedef typename Rep::pointer               pointer;    //!< pointer type to elements
		typedef const pointer                 const_pointer;    //!< const pointer type
		typedef BlasMatrix<_Field>                   Self_t;    //!< Self type

	protected:
		size_t			_row;
		size_t			_col;
		Rep			_rep;
		pointer			_ptr;
		const Field		& _F;
		MatrixDomain<Field>      _MD;
		VectorDomain<Field>      _VD;
		bool		 _use_fflas ;


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
		void createBlasMatrix (const BlasMatrix<_Field> & A) ;

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

		/*!@internal constructor from vector of elements.
		 * @param v pointer to \c Element s
		 */
		void createBlasMatrix ( const Element * v) ;

		/*!@internal constructor from vector of elements.
		 * @param v std::vector of \c Element s
		 */
		void createBlasMatrix ( const std::vector<Element> & v) ;
		/*! @internal
		 * @}
		 */

	public:

		//////////////////
		// CONSTRUCTORS //
		//////////////////


		/*! Allocates a new \f$ 0 \times 0\f$ matrix.
		*/
		BlasMatrix (const _Field &F) ;

		/*! Allocates a new \f$ 0 \times 0\f$ matrix.
		*/
		BlasMatrix () ;

		/*! Allocates a new \f$ m \times n\f$ matrix.
		 * @param m rows
		 * @param n cols
		 */
		//@{

		template<class T>
		BlasMatrix (const _Field &F, uint64_t m, T n) ;

		template<class T>
		BlasMatrix (const _Field &F, int64_t m, T n) ;

		template<class T>
		BlasMatrix (const _Field &F, uint32_t m, T  n) ;

		template<class T>
		BlasMatrix (const _Field &F, int32_t m, T n) ;

		template<class T>
		BlasMatrix (const _Field &F, Integer & m, T n) ;

		//@}


		/*! Constructor from a matrix stream.
		 * @param ms matrix stream.
		 */
		BlasMatrix(MatrixStream<_Field>& ms) ;

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
		template<class _Matrix>
		BlasMatrix (const _Matrix &A,  const _Field &F) ;

		/*! Copy Constructor of a matrix (copying data).
		 * @param A matrix to be copied.
		 */
		BlasMatrix (const BlasMatrix<_Field>& A) ;

		/*- Copy Constructor of a matrix (copying data).
		 * @param A matrix to be copied.
		 */
		// BlasMatrix (const BlasSubmatrix<_Field>& A) ;

		/*! Create a BlasMatrix from a vector of elements
		 * @param A matrix to be copied.
		 */
		BlasMatrix (const _Field &F, const std::vector<Element>& v,
			    size_t m, size_t n) ;

		/*! Create a BlasMatrix from an array of elements
		 * @param A matrix to be copied.
		 */
		BlasMatrix (const _Field &F, const Element * v,
			    size_t m, size_t n) ;


		/** Constructor using a finite vector stream (stream of the rows).
		 * @param  F The field of entries; passed so that arithmetic may be done
		 *           on elements
		 * @param  stream A vector stream to use as a source of vectors for this
		 *                matrix
		 */
		template <class StreamVector>
		BlasMatrix (const Field &F, VectorStream<StreamVector> &stream) ;

		/// Destructor.
		~BlasMatrix () ;

		//! operator = (copying data)
		BlasMatrix<_Field>& operator= (const BlasMatrix<_Field>& A) ;

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
		BlasMatrix<_Field> transpose(BlasMatrix<_Field> & tM) const ;


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
		// template <class Field>
		std::istream &read (std::istream &file/*  , const Field &F*/);

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

		///////////////////
		//   BLACK BOX   //
		///////////////////


		template <class Vector1, class Vector2>
		Vector1&  apply (Vector1& y, const Vector2& x) const ;

		template <class Vector1, class Vector2>
		Vector1&  applyTranspose (Vector1& y, const Vector2& x) const ;

		const _Field& field() const;
		_Field& field() ;
		void setField(const _Field & F) { _F = F ; };


	}; // end of class BlasMatrix

	template <class _Field>
	struct MatrixTraits< BlasMatrix<_Field> > {
		typedef BlasMatrix<_Field> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};

	template <class _Field>
	struct MatrixTraits< const BlasMatrix<_Field> > {
		typedef const BlasMatrix<_Field> MatrixType;
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
	 * A @ref BlasSubmatrix is a matrix of \p _Field::Element, with the structure of BLAS matrices.
	 * It is basically a read/write view on a vector of \p _Field::Element.
	 * In the Mother model, a @ref BlasSubmatrix is not allocated.
	 * <p>
	 * This matrix type conforms to the same interface as @ref BlasMatrix,
	 * except that you cannot resize it. It represents a submatrix of a dense
	 * matrix. Upon construction, one can freely manipulate the entries in the
	 * DenseSubmatrix, and the corresponding entries in the underlying
	 * @ref BlasMatrix will be modified.


	 */
	template <class _Field>
	class BlasSubmatrix {
	public :
		typedef _Field                       Field;
		typedef typename Field::Element      Element;      //!< Element type
		typedef BlasSubmatrix<_Field>        Self_t;       //!< Self type


	protected:
		BlasMatrix<_Field> *_M;       //!< Parent BlasMatrix (ie raw vector)
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
		BlasSubmatrix (const BlasMatrix<_Field> &M,
			       size_t rowbeg,
				size_t colbeg,
				size_t Rowdim,
				size_t Coldim);

		/** Constructor from an existing @ref DenseMatrixBase
		 * \param M Pointer to @ref DenseMatrixBase of which to construct submatrix
		 */
		BlasSubmatrix (const BlasMatrix<_Field> &M);


		/** Constructor from an existing submatrix and dimensions
		 * @param SM Constant reference to BlasSubmatrix from which to
		 *           construct submatrix
		 * @param rowbeg Starting row
		 * @param colbeg Starting column
		 * @param Rowdim Row dimension
		 * @param Coldim Column dimension
		 */
		BlasSubmatrix (const BlasSubmatrix<_Field> &SM,
				size_t rowbeg,
				size_t colbeg,
				size_t Rowdim,
				size_t Coldim);

		/** Copy constructor.
		 * @param SM Submatrix to copy
		 */
		BlasSubmatrix (const BlasSubmatrix<_Field> &SM);


		/*  Members  */

		/** Assignment operator.
		 * Assign the given submatrix to this one
		 * This is <i>only</i> renaming !
		 * There is no copy because BlasSubmatrix owns nothing.
		 * @param SM Submatrix to assign
		 * @return Reference to this submatrix
		 */
		BlasSubmatrix &operator = (const BlasSubmatrix<_Field> &SM);

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
		// template<class Field>
		std::istream& read (std::istream &file/*, const Field& field*/);


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
		typedef typename BlasMatrix<_Field>::RowIterator            RowIterator;
		typedef typename BlasMatrix<_Field>::ConstRowIterator       ConstRowIterator;
		typedef typename BlasMatrix<_Field>::Row                    Row;
		typedef typename BlasMatrix<_Field>::ConstRow               ConstRow;
		//@} Row Iterators

		/** @name typedef'd Column Iterators.
		 *\brief
		 * The columns iterator gives the columns of the
		 * matrix in ascending order. Dereferencing the iterator yields
		 * a column vector in dense format
		 * @{
		 */
		typedef typename BlasMatrix<_Field>::ColIterator            ColIterator;
		typedef typename BlasMatrix<_Field>::ConstColIterator       ConstColIterator;
		typedef typename BlasMatrix<_Field>::Col                    Col;
		typedef typename BlasMatrix<_Field>::Column                 Column;
		typedef typename BlasMatrix<_Field>::ConstCol               ConstCol;
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

		///////////////////
		//   BLACK BOX   //
		///////////////////


		template <class Vector1, class Vector2>
		Vector1&  apply (Vector1& y, const Vector2& x) const
		{
			//_stride ?
			if (_M->_use_fflas){
				//!@bug this supposes &x[0]++ == &x[1]
				FFLAS::fgemv( _M->_F, FFLAS::FflasNoTrans,
					      _row, _col,
					      _M->_F.one,
					      _M->_ptr, getStride(),
					      &x[0],1,
					      _M->_F.zero,
					      &y[0],1);
			}
			else {
				_M->_MD. vectorMul (y, *this, x);
#if 0
				typename BlasMatrix<_Field>::ConstRowIterator i = this->rowBegin ();
				typename Vector1::iterator j = y.begin ();

				for (; j != y.end (); ++j, ++i)
					_VD.dot (*j, *i, x);
#endif
			}
			return y;
		}

		template <class Vector1, class Vector2>
		Vector1&  applyTranspose (Vector1& y, const Vector2& x) const
		{

			//_stride ?
			if (_M->_use_fflas) {
				FFLAS::fgemv( _M->_F, FFLAS::FflasTrans,
					      _row, _col,
					      _M->_F.one,
					      _M->_ptr, getStride(),
					      &x[0],1,
					      _M->_F.zero,
					      &y[0],1);
			}
			else {
				typename BlasMatrix<_Field>::ConstColIterator i = this->colBegin ();
				typename Vector1::iterator j = y.begin ();
				for (; j != y.end (); ++j, ++i)
					_M->_VD.dot (*j, x, *i);
			}

			return y;
		}

	};

	template <class _Field>
	struct MatrixTraits< BlasSubmatrix<_Field> > {
		typedef BlasSubmatrix<_Field> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};

	template <class _Field>
	struct MatrixTraits< const BlasSubmatrix<_Field> > {
		typedef const BlasSubmatrix<_Field> MatrixType;
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
	template <class _Field>
	class TriangularBlasMatrix: public BlasMatrix<_Field> {

	protected:

		LinBoxTag::Shape          _uplo; //!< upper or lower triangular
		LinBoxTag::Diag           _diag; //!< unit or non unit diagonal

	public:
		typedef _Field                       Field;
		typedef typename Field::Element      Element;      //!< Element type
		typedef BlasMatrix<_Field>           Father_t;
		typedef TriangularBlasMatrix<_Field> Self_t;


		/*! Constructor for a new \c TriangularBlasMatrix.
		 * @param m rows
		 * @param n cols
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		TriangularBlasMatrix (const Field & F,
				      const size_t m, const size_t n,
				      LinBoxTag::Shape x=LinBoxTag::Upper,
				      LinBoxTag::Diag y= LinBoxTag::NonUnit) ;

		/*! Constructor from a \c BlasMatrix (copy).
		 * @param A matrix
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		TriangularBlasMatrix (const BlasMatrix<_Field>& A,
				      LinBoxTag::Shape x=LinBoxTag::Upper,
				      LinBoxTag::Diag y= LinBoxTag::NonUnit) ;

		/*! Constructor from a \c BlasMatrix (no copy).
		 * @param A matrix
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		TriangularBlasMatrix (BlasMatrix<_Field>& A,
				      LinBoxTag::Shape x=LinBoxTag::Upper,
				      LinBoxTag::Diag y= LinBoxTag::NonUnit) ;

		/*! Constructor from a \c TriangularBlasMatrix (copy).
		 * @param A matrix
		 */
		TriangularBlasMatrix (const TriangularBlasMatrix<_Field>& A) ;

		/*! Generic constructor from a \c Matrix (no copy).
		 * @param A matrix
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		template<class Matrix>
		TriangularBlasMatrix (const Matrix& A,
				      LinBoxTag::Shape x=LinBoxTag::Upper,
				      LinBoxTag::Diag y= LinBoxTag::NonUnit) ;

		/// get the shape of the matrix (upper or lower)
		LinBoxTag::Shape getUpLo() const ;

		/// Is the diagonal implicitly unit ?
		LinBoxTag::Diag getDiag() const ;

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

#include "linbox/field/multimod-field.h"
namespace LinBox
{
	template<>
	class BlasMatrix<MultiModDouble> {

	public:

		typedef MultiModDouble         Field;
		typedef std::vector<double>  Element;
		typedef BlasMatrix<MultiModDouble> Self_t;

	protected:

		MultiModDouble                 _F;
		const std::vector<MatrixDomain<Modular<double> > >   _MD;
		size_t                  _row,_col;
		Element                _One,_Zero;
		std::vector<BlasMatrix<Modular<double> >* > _rep;
		std::vector<double>       _entry;
	public:


		//BlasMatrix () {}

		BlasMatrix (const MultiModDouble& F) :
			_F(F) , _rep(F.size()), _entry(F.size())
		{}

		BlasMatrix (const Field& F, size_t m, size_t n, bool alloc=true) :
			_F(F), _row(m) , _col(n) , _rep(F.size()),  _entry(F.size())
		{
			for (size_t i=0;i<_rep.size();++i)
				_rep[i] =  new BlasMatrix<Modular<double> > (F.getBase(i), m, n);
		}

		BlasMatrix (const BlasMatrix<MultiModDouble> & A):
			_F(A._F),_row(A._row), _col(A._col),
			_rep(A._rep.size()), _entry(A._entry)
		{

			for (size_t i=0;i<_rep.size();++i)
				_rep[i]= new  BlasMatrix<Modular<double> > (const_cast<BlasMatrix<Modular<double> >& >( *A._rep[i]));
		}


		const BlasMatrix<MultiModDouble>& operator=(const BlasMatrix<MultiModDouble> & A)
		{
			_F   = A._F;
			_row = A._row;
			_col = A._col;
			_rep = std::vector<BlasMatrix<Modular<double> >* >(A._rep.size());
			_entry = A._entry;
			for (size_t i=0;i<_rep.size();++i)
				_rep[i]= new  BlasMatrix<Modular<double> > (const_cast<BlasMatrix<Modular<double> >& >( *A._rep[i]));
			return *this;
		}


		~BlasMatrix() {for (size_t i=0; i< _rep.size();++i) {delete _rep[i];} }

		template <class Vector1, class Vector2>
		Vector1&  apply (Vector1& y, const Vector2& x) const
		{
			for (size_t i=0;i<_rep.size();++i) {
				std::vector<double> x_tmp(x.size()), y_tmp(y.size());
				for (size_t j=0;j<x.size();++j)
					x_tmp[j]= x[j][i];

				_rep[i]->apply(y_tmp, x_tmp);

				for (size_t j=0;j<y.size();++j){
					y[j][i]=y_tmp[j];

				}
			}

			return y;
		}

		template <class Vector1, class Vector2>
		Vector1&  applyTranspose (Vector1& y, const Vector2& x) const
		{
			for (size_t i=0;i<_rep.size();++i) {
				std::vector<double> x_tmp(x.size()), y_tmp(y.size());
				for (size_t j=0;j<x.size();++j)
					x_tmp[i]= x[j][i];

				_rep[i]->applyTranspose(y_tmp, x_tmp);

				for (size_t j=0;j<y.size();++j)
					y[j][i]=y_tmp[i];
			}

			return y;
		}

#if 0
		template<typename _Tp1>
		struct rebind
		{
			typedef BlasMatrix<_Tp1> other;

			void operator() (other *& Ap, const Self_t& A, const _Tp1& F) {
				Ap = new other(F, A.rowdim(), A.coldim());
				Hom<Field, _Tp1> hom(A. field(), F);

				hom.image (*Ap_p, *A_p);
			}
		};
#endif

		size_t rowdim() const {return _row;}

		size_t coldim() const {return _col;}


		const Field &field() const  {return _F;}


		std::ostream& write(std::ostream& os) const
		{
			for (size_t i=0;i<_rep.size();++i)
				_rep[i]->write(os);
			return os;
		}


		void setEntry (size_t , size_t j, const Element &a_ij)
		{
			for (size_t i=0; i< _rep.size();++i)
				_rep[i]->setEntry(i,j,a_ij[i]);
		}


		const Element& getEntry (size_t , size_t j)
		{
			for (size_t i=0; i< _rep.size();++i)
				_entry[i]=_rep[i]->getEntry(i,j);
			return _entry;
		}

		BlasMatrix<Modular<double> >*& getMatrix(size_t i) {return _rep[i];}

	};

}

#include "blas-matrix.inl"
#include "blas-submatrix.inl"
#include "blas-triangularmatrix.inl"

#endif // __LINBOX_blas_matrix_H

