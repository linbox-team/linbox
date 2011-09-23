/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/matrix/dense.h
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
 * Move from blackbox/dense-base.h to matrix/dense.h
 * --------------------------------------------------------
 * 2002-11-29  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Swap the order of arguments in read and write, so that it is consistent with
 * SparseMatrixBase
 * --------------------------------------------------------
 * 2002-10-28  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Rename ColOfRowsIterator as RowIterator; similarly with RowOfColsIterator
 * --------------------------------------------------------
 * 2002-10-27  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Split out container/iterator functionality into DenseMatrixBase
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

/*! @file matrix/dense.h
 * @ingroup matrix
 * @brief Blackbox dense matrix.
 */

#ifndef __LINBOX_matrix_dense_H
#define __LINBOX_matrix_dense_H

#include <iostream>
#include <vector>
#include <fstream>

#include "linbox/vector/subiterator.h"
#include "linbox/vector/subvector.h"
#include "linbox/vector/stream.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/linbox-config.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/field/hom.h"


namespace LinBox
{

	namespace Protected
	{
		/** @brief LinBox dense matrix template.
		 * This is a class of dense matrices
		 * templatized by the entry type, the \c Element type of some @link Fields field@endlink.
		 * The matrix is stored as a one dimensional STL \c std::vector of the elements, by rows.
		 * The interface provides for iteration over rows and over columns.
		 *
		 * The class \c LinBox::DenseMatrix builds on this base.
		 *
		 * Currently, only dense vectors are supported when doing matrix-vector applies.
		 *
		 * \ingroup matrix
		 */
		template <class _Element>
		class DenseMatrixBase {
		public:

			typedef _Element                            Element; //!< Element type
			typedef typename RawVector<Element>::Dense      Rep; //!< No doc
			typedef DenseMatrixBase<_Element>            Self_t; //!< Self type

			template<typename _Tp1>
			struct rebind {
				typedef DenseMatrixBase<typename _Tp1::Element> other;
			};

			/// NULL constructor.
			DenseMatrixBase () :
				_rows (0), _cols (0)
			{}

			/** Constructor.
			 * @param  m  row dimension
			 * @param  n  column dimension
			 */
			DenseMatrixBase (size_t m, size_t n) :
				_rep (m * n), _rows (m), _cols (n)
			{}

			/** Constructor.
			 * @param data vector to take data from (such that <code>A[i,j] == data[i*n+j]</code>)
			 * @param  m  row dimension
			 * @param  n  column dimension
			 */
			DenseMatrixBase (std::vector<Element> & data,
					 size_t m, size_t n) :
				_rep (data), _rows (m), _cols (n)
			{}

			/** Constructor from a matrix stream.
			 * @param ms matrix stream
			 */
			template< class Field >
			DenseMatrixBase( MatrixStream<Field>& ms );

			/*! copy constructor.
			 * @param M Dense matrix to be copied.
			 */
			DenseMatrixBase (const DenseMatrixBase &Mat) :
				_rep (Mat._rep),_rows (Mat._rows), _cols (Mat._cols)
			{}

			//! destructor.
			//~DenseMatrixBase(){}

			/*! copy operator.
			 * @param M matrix to be copied.
			 */
			DenseMatrixBase& operator= (const DenseMatrixBase& Mat)
			{
				(*this)._rep  = Mat._rep;
				(*this)._rows = Mat._rows;
				(*this)._cols = Mat._cols;
				return (*this);
			}

			/** Get a pointer on the storage of the elements.
			 * @returns a pointer on Elements
			 /todo What is this?
			 */
			Element* FullIterator() const
			{
				return const_cast<Element*>(&_rep[0]);
			}

			/** Get the number of rows in the matrix
			 * @returns Number of rows in matrix
			 */
			size_t rowdim () const
			{
				return _rows;
			}

			/** Get the number of columns in the matrix.
			 * @returns Number of columns in matrix
			 */
			size_t coldim () const
			{
				return _cols;
			}

			/** Resize the matrix to the given dimensions.
			 * The state of the matrix's entries after a call to this method is
			 * undefined
			 * @param m Number of rows
			 * @param n Number of columns
			 * @param val
			 */
			void resize (size_t m, size_t n, const Element& val = Element())
			{
				_rows = m;
				_cols = n;
				_rep.resize (m * n, val);
			}

			/** Read the matrix from an input stream.
			 * The stream is in SMS or DENSE format
			 * @param file Input stream from which to read
			 * @param F Field over which to read
			 */
			template <class Field>
			std::istream &read (std::istream &file, const Field &F);

			/** Write the matrix to an output stream.
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


			/** Set the entry at the (i, j) position to a_ij.
			 * @param i Row number, 0...rowdim () - 1
			 * @param j Column number 0...coldim () - 1
			 * @param a_ij Element to set
			 */
			void setEntry (size_t i, size_t j, const Element &a_ij)
			{
				_rep[i * _cols + j] = a_ij;
			}

			/** Get a writeable reference to the entry in the (i, j) position.
			 * @param i Row index of entry
			 * @param j Column index of entry
			 * @returns Reference to matrix entry
			 */
			Element &refEntry (size_t i, size_t j)
			{
				return _rep[i * _cols + j];
			}

			/** Get a read-only reference to the entry in the (i, j) position.
			 * @param i Row index
			 * @param j Column index
			 * @returns Const reference to matrix entry
			 */
			const Element &getEntry (size_t i, size_t j) const
			{
				return _rep[i * _cols + j];
			}

			/** Copy the (i, j) entry into x, and return a reference to x.
			 * This form is more in the Linbox style and is provided for interface
			 * compatibility with other parts of the library
			 * @param x Element in which to store result
			 * @param i Row index
			 * @param j Column index
			 * @returns Reference to x
			 */
			Element &getEntry (Element &x, size_t i, size_t j) const
			{
				x = _rep[i * _cols + j]; return x;
			}


			/** @name Column of rows iterator
			 * \brief
			 * The column of rows iterator traverses the rows of the
			 * matrix in ascending order. Dereferencing the iterator yields
			 * a row vector in dense format
			 */
			//@{
			typedef Subvector<typename Rep::iterator, typename Rep::const_iterator> Row;
			typedef Subvector<typename Rep::const_iterator>                    ConstRow;

			class RowIterator;
			class ConstRowIterator;

			RowIterator rowBegin ();
			RowIterator rowEnd ();
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

			class ColIterator;
			class ConstColIterator;

			ColIterator colBegin ();
			ColIterator colEnd ();
			ConstColIterator colBegin () const;
			ConstColIterator colEnd () const;
			//@}

			/** @name Raw iterator
			 * \brief
			 *
			 * The raw iterator is a method for accessing all entries in the matrix
			 * in some unspecified order. This can be used, e.g. to reduce all
			 * matrix entries modulo a prime before passing the matrix into an
			 * algorithm.
			 */
			//@{
			typedef typename Rep::iterator RawIterator;
			typedef typename Rep::const_iterator ConstRawIterator;

			RawIterator rawBegin ();
			RawIterator rawEnd   ();
			ConstRawIterator rawBegin () const;
			ConstRawIterator rawEnd   () const;
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
			class RawIndexedIterator;
			class ConstRawIndexedIterator;

			RawIndexedIterator rawIndexedBegin ();
			RawIndexedIterator rawIndexedEnd   ();
			ConstRawIndexedIterator rawIndexedBegin () const;
			ConstRawIndexedIterator rawIndexedEnd   () const;
			//@}

			/** Retrieve a reference to a row.
			 * Since rows may also be indexed, this allows A[i][j] notation
			 * to be used.
			 * @param i Row index
			 */
			//@{
			Row operator[] (size_t i)
			{
				return Row (_rep.begin () + i * _cols, _rep.begin () + i * _cols + _cols);
			}

			ConstRow operator[] (size_t i) const
			{
				return Row (_rep.begin () + i * _cols, _rep.begin () + i * _cols + _cols);
			}
			//@}

			/** Compute column density.
			 * @param v
			 */
			template <class Vector>
			Vector &columnDensity (Vector &v) const
			{
				std::fill (v.begin (), v.end (), _rows);
				return v;
			}

		protected:
			std::vector<Element>  _rep;
			size_t                _rows, _cols;

		}; //class DenseMatrixBase

	} // Protected

	/*! Write a matrix to a stream.
	 * The C++ way using <code>operator<<</code>
	 * @param o output stream
	 * @param M matrix to write.
	 */
	template<class T>
	std::ostream& operator<< (std::ostream & o, const Protected::DenseMatrixBase<T> & Mat)
	{
		return Mat.write(o);
	}

	template <class Element>
	struct MatrixTraits< Protected::DenseMatrixBase<Element> > {
		typedef Protected::DenseMatrixBase<Element> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};

	template <class Element>
	struct MatrixTraits< const Protected::DenseMatrixBase<Element> > {
		typedef const Protected::DenseMatrixBase<Element> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};

} // namespace LinBox

#include "dense.inl"

#endif // __LINBOX_matrix_dense_H

