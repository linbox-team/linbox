/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/matrix/dense-submatrix.h
 * Copyright (C) 2001 B. David Saunders,
 *               2001-2002 Bradford Hovinen,
 *               2002 Zhendong Wan
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>,
 *            Zhendong Wan <wan@mail.eecis.udel.edu>
 *
 * evolved from dense-matrix.h by -bds, Zhendong Wan
 *
 * --------------------------------------------------------
 * 2003-01-11  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Move from blackbox to matrix
 * -----------------------------------------------------------
 * 2002-11-30  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Have the constructor take a reference rather than a pointer
 * -----------------------------------------------------------
 * 2002-10-27  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Rename from densesubmatrix.h
 *
 * Constructor modifications: changed the interface to match Submatrix
 *
 * Don't parameterize by Field, but instead by Element; remove all the black box
 * apply stuff
 * -----------------------------------------------------------
 *
 * See COPYING for license information
 */

/*! @file matrix/dense-submatrix.h
 * @ingroup matrix
 * @brief Representation of a submatrix of a dense matrix, not resizeable.
 * This matrix type conforms to the \c LinBox::DenseMatrixBase interface.
 * \c LinBox::BlasMatrix is an example of DenseSubmatrix.
 */

#ifndef __LINBOX_dense_submatrix_H
#define __LINBOX_dense_submatrix_H

#include "linbox/linbox-config.h"

#include "linbox/util/debug.h"
#include "linbox/matrix/dense.h"
#include "linbox/matrix/matrix-domain.h"

namespace LinBox
{

	namespace Protected
	{
	/** @brief %Submatrix of a dense matrix
	 *
	 * This matrix type conforms to the same interface as @ref DenseMatrixBase,
	 * except that you cannot resize it. It represents a submatrix of a dense
	 * matrix. Upon construction, one can freely manipulate the entries in the
	 * DenseSubmatrix, and the corresponding entries in the underlying
	 * DenseMatrixBase will be modified.

	 \ingroup matrix
	 */
	template<class _Element>
	class DenseSubmatrix {
	public:

		/*  Iterators */

		//! @name Forward declaration of Raw Iterators.
		//@{
		// RawIterator.
		class RawIterator  ;
		// Const RawIterator.
		class ConstRawIterator ;

		 // RawIndexed
		class RawIndexedIterator ;
		// Const RawIndexed
		class ConstRawIndexedIterator ;
		//@}

		typedef _Element                  Element;       //!< Element type
		typedef DenseSubmatrix<_Element>   Self_t;       //!< Self type


		/** @name typedef'd Row Iterators.
		 *\brief
		 * The row iterator gives the rows of the
		 * matrix in ascending order. Dereferencing the iterator yields
		 * a row vector in dense format
		 * @{
		 */
		typedef typename DenseMatrixBase<Element>::RowIterator            RowIterator;
		typedef typename DenseMatrixBase<Element>::ConstRowIterator       ConstRowIterator;
		typedef typename DenseMatrixBase<Element>::Row                    Row;
		typedef typename DenseMatrixBase<Element>::ConstRow               ConstRow;
		 //@} Row Iterators

		/** @name typedef'd Column Iterators.
		 *\brief
		 * The columns iterator gives the columns of the
		 * matrix in ascending order. Dereferencing the iterator yields
		 * a column vector in dense format
		 * @{
		 */
		typedef typename DenseMatrixBase<Element>::ColIterator            ColIterator;
		typedef typename DenseMatrixBase<Element>::ConstColIterator       ConstColIterator;
		typedef typename DenseMatrixBase<Element>::Col                    Col;
		typedef typename DenseMatrixBase<Element>::Column                 Column;
		typedef typename DenseMatrixBase<Element>::ConstCol               ConstCol;
		//@} // Column Iterators


		/*  constructors */

		/** NULL constructor.  */
		DenseSubmatrix () :
			_M(NULL)
	       	{}

		/** Constructor from an existing @ref DenseMatrixBase  and dimensions.
		 * \param M Pointer to @ref DenseMatrixBase of which to construct submatrix
		 * \param row Starting row
		 * \param col Starting column
		 * \param rowdim Row dimension
		 * \param coldim Column dimension
		 */
		DenseSubmatrix (DenseMatrixBase<Element> &M,
				size_t row,
				size_t col,
				size_t rowdim,
				size_t coldim);

		/** Constructor from an existing @ref DenseMatrixBase
		 * \param M Pointer to @ref DenseMatrixBase of which to construct submatrix
		 */
		DenseSubmatrix (DenseMatrixBase<Element> &M);


		/** Constructor from an existing submatrix and dimensions
		 * @param SM Constant reference to DenseSubmatrix from which to
		 *           construct submatrix
		 * @param row Starting row
		 * @param col Starting column
		 * @param rowdim Row dimension
		 * @param coldim Column dimension
		 */
		DenseSubmatrix (const DenseSubmatrix<Element> &SM,
				size_t row,
				size_t col,
				size_t rowdim,
				size_t coldim);

		/** Copy constructor.
		 * @param SM Submatrix to copy
		 */
		DenseSubmatrix (const DenseSubmatrix<Element> &SM);

		/*  Members  */

		/** Assignment operator.
		 * Assign the given submatrix to this one
		 * @param SM Submatrix to assign
		 * @return Reference to this submatrix
		 */
		DenseSubmatrix &operator = (const DenseSubmatrix<Element> &SM);

		/** Get the number of rows in the matrix
		 * @return Number of rows in matrix
		 */
		size_t rowdim () const
		{
			return _end_row - _beg_row;
		}

		/** Get the number of columns in the matrix
		 * @return Number of columns in matrix
		 */
		size_t coldim () const
		{
			return _end_col - _beg_col;
		}

		template<typename _Tp1>
		struct rebind {
			typedef DenseSubmatrix<typename _Tp1::Element> other;
		};


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
				     bool mapleFormat = false) const;

		/** Write the matrix to an output stream.
		 * This a raw version of \c write(os,F) (no field is given).
		 * @param os Output stream to which to write
		 * @param mapleFormat write in Maple(r) format ?
		 */
		std::ostream& write (std::ostream &os,
				     bool mapleFormat = false) const;


		/** Set the entry at (i, j).
		 * @param i Row number, 0...rowdim () - 1
		 * @param j Column number 0...coldim () - 1
		 * @param a_ij Element to set
		 */
		void setEntry (size_t i, size_t j, const Element &a_ij)
		{
			_M->setEntry (_beg_row + i, _beg_col + j, a_ij);
		}

		/** Get a writeable reference to an entry in the matrix.
		 * @param i Row index of entry
		 * @param j Column index of entry
		 * @return Reference to matrix entry
		 */
		Element &refEntry (size_t i, size_t j)
		{
			return _M->refEntry (i + _beg_row, j + _beg_col);
		}

		/** Get a read-only individual entry from the matrix.
		 * @param i Row index
		 * @param j Column index
		 * @return Const reference to matrix entry
		 */
		const Element &getEntry (size_t i, size_t j) const
		{
			return _M->getEntry (i + _beg_row, j + _beg_col);
		}

		/** Get an entry and store it in the given value.
		 * This form is more in the Linbox style and is provided for interface
		 * compatibility with other parts of the library
		 * @param x Element in which to store result
		 * @param i Row index
		 * @param j Column index
		 * @return Reference to x
		 */
		Element &getEntry (Element &x, size_t i, size_t j) const
		{
			return _M->getEntry (x, i + _beg_row, j + _beg_col);
		}

#if 0 /*  craquage */
		Element & operator[] (size_t i, size_t j)
		{
			return refEntry(i,j);
		}

		const Element & operator[] (size_t i, size_t j) const
		{
			return getEntry(i,j);
		}
#endif

		/// iterator to the begining of a row
		RowIterator rowBegin ();
		/// iterator to the end of a row
		RowIterator rowEnd ();
		/// const iterator to the begining of a row
		ConstRowIterator rowBegin () const;
		/// const iterator to the end of a row
		ConstRowIterator rowEnd () const;

		ColIterator colBegin ();
		ColIterator colEnd ();
		ConstColIterator colBegin () const;
		ConstColIterator colEnd () const;

		RawIterator rawBegin ();
		RawIterator rawEnd ();
		ConstRawIterator rawBegin () const;
		ConstRawIterator rawEnd () const;


		RawIndexedIterator rawIndexedBegin();
		RawIndexedIterator rawIndexedEnd();
		ConstRawIndexedIterator rawIndexedBegin() const;
		ConstRawIndexedIterator rawIndexedEnd() const;

#if 0 /*  operator[] */
		/*- Retrieve a reference to a row
		 * @param i Row index
		 */
		Row operator[] (int i);               not actually used, causes a compile error...
		ConstRow operator[] (int i) const;
#endif

		/*! Creates a transposed matrix of \c *this.
		 * @param[in] tM
		 * @return the transposed matrix of this.
		 */
		DenseSubmatrix<Element> transpose(DenseMatrixBase<Element> & tM)
		{
			linbox_check(tM.coldim() == rowdim());
			linbox_check(tM.rowdim() == coldim());
			// DenseMatrixBase<Element> tM(coldim(),rowdim());
			DenseSubmatrix<Element>  tA(tM);
			for (size_t i = 0 ; i < rowdim(); ++i)
				for (size_t j = 0 ; j < coldim(); ++j)
					tA.setEntry(j,i,getEntry(i,j));
			return tA;
		}

		/*! Creates a transposed matrix of \c *this.
		 * @return the transposed matrix of this.
		 */
		DenseSubmatrix<Element> & transpose(DenseSubmatrix<Element> & tA)
		{
			for (size_t i = 0 ; i < rowdim(); ++i)
				for (size_t j = 0 ; j < coldim(); ++j)
					tA.setEntry(j,i,getEntry(i,j));
			return tA;
		}

	protected:
		DenseMatrixBase<Element> *_M;
		size_t _beg_row;
		size_t _end_row;
		size_t _beg_col;
		size_t _end_col;
	};
	} // Protected

	/*! Write a matrix to a stream.
	 * The C++ way using <code>operator<<</code>
	 * @param o output stream
	 * @param Mat matrix to write.
	 */
	template<class T>
	std::ostream& operator<< (std::ostream & o, const Protected::DenseSubmatrix<T> & Mat)
	{
		return Mat.write(o);
	}


	/*! @internal
	 * @brief MatrixTraits
	 */
	template <class Element>
	struct MatrixTraits< Protected::DenseSubmatrix<Element> > {
		typedef Protected::DenseSubmatrix<Element> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};

} // namespace LinBox

#include "linbox/matrix/dense-submatrix.inl"

#endif // __LINBOX_dense_submatrix_H


