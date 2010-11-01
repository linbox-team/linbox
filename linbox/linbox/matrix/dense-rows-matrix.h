/* linbox/matrix/dense-rows-matrix.h
 * Copyright (C) 2003 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 * 
 * --------------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_dense_rows_matrix_H
#define __LINBOX_dense_rows_matrix_H

#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>

#include "linbox/blackbox/factory.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/matrix/matrix-domain.h"

namespace LinBox
{

/** Dense row-wise matrix container
 *
 * This class implements a dense matrix, storing the data as a vector of vectors
 * of the given type, in the same manner as @ref{SparseMatrixBase}. It provides
 * only row iterators.
 *
 * @param Vector Row vector type
 */
template <class _Row>
class DenseRowsMatrix
{
    public:

	typedef _Row Row;
	typedef typename Row::value_type Element;
	typedef typename std::vector<Row> Rep;

	/** Constructor.
	 * Note: the copy constructor and operator= will work as intended
	 *       because of STL's container design
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	DenseRowsMatrix (size_t m, size_t n);

	/** Copy constructor.
	 */
	DenseRowsMatrix (const DenseRowsMatrix &A)
		: _A (A._A), _m (A._m), _n (A._n) 
	{}

	/** Destructor. */
	~DenseRowsMatrix () {}

	/** Retreive row dimension of the matrix.
	 * @return integer number of rows of DenseRowsMatrix matrix.
	 */
	size_t rowdim () const { return _m; }

	/** Retreive column dimension of matrix.
	 * @return integer number of columns of DenseRowsMatrix matrix.
	 */
	size_t coldim () const { return _n; }

	/** @name Input and output
	 */
	//@{

	/** Read a matrix from the given input stream using field read/write
	 * @param is Input stream from which to read the matrix
	 * @param F Field with which to read
	 * @param format Format of input matrix
	 */
	template <class Field>
	std::istream &read (std::istream &is, const Field &F);

	/** Read a matrix from the given input stream using standard operators
	 * @param is Input stream from which to read the matrix
	 * @param format Format of input matrix
	 */
	std::istream &read (std::istream &is);

	/** Write a matrix to the given output stream using field read/write
	 * @param os Output stream to which to write the matrix
	 * @param F Field with which to write
	 * @param format Format with which to write
	 */
	template <class Field>
	std::ostream &write (std::ostream &os, const Field &F) const;

	/** Write a matrix to the given output stream using standard operators
	 * @param os Output stream to which to write the matrix
	 * @param format Format with which to write
	 */
	std::ostream &write (std::ostream &os) const;

	//@}

	/** @name Access to matrix elements
	 */
	//@{

	/** Set an individual entry
	 * Setting the entry to 0 will remove it from the matrix
	 * @param i Row index of entry
	 * @param j Column index of entry
	 * @value Value of the new entry
	 */
	void setEntry (size_t i, size_t j, const Element &value)
		{ _A[i][j] = value; }

	/** Get a writeable reference to an entry in the matrix
	 * If there is no entry at the position (i, j), then a new entry
	 * with a value of zero is inserted and a reference  to it is
	 * returned.
	 * @param i Row index of entry
	 * @param j Column index of entry
	 * @return Reference to matrix entry
	 */
	Element &refEntry (size_t i, size_t j)
		{ return _A[i][j]; }

	/** Get a read-only individual entry from the matrix
	 * @param i Row index
	 * @param j Column index
	 * @return Const reference to matrix entry
	 */
	const Element &getEntry (size_t i, size_t j) const
		{ return _A[i][j]; }

	/** Get an entry and store it in the given value
	 * This form is more in the Linbox style and is provided for interface
	 * compatibility with other parts of the library
	 * @param x Element in which to store result
	 * @param i Row index
	 * @param j Column index
	 * @return Reference to x
	 */
	Element &getEntry (Element &x, size_t i, size_t j) const
		{ return x = _A[i][j]; }

	/** @name Columns of rows iterator
	 * The columns of row iterator gives each of the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in sparse sequence format
	 */

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	RowIterator rowBegin () { return _A.begin (); }
	RowIterator rowEnd () { return _A.end (); }
	ConstRowIterator rowBegin () const { return _A.begin (); }
	ConstRowIterator rowEnd () const { return _A.end (); }

	/** @name Raw iterator
	 * The raw iterator is a method for accessing all nonzero
	 * entries in the matrix in some unspecified order. This can be
	 * used, e.g. to reduce all matrix entries modulo a prime before
	 * passing the matrix into an algorithm.
	 */

	class RawIterator;
	class ConstRawIterator;

	RawIterator rawBegin ();
	RawIterator rawEnd ();
	ConstRawIterator rawBegin () const;
	ConstRawIterator rawEnd () const;

	/** @name Index iterator
	 * The index iterator gives the row, column indices of all matrix
	 * elements in the same order as the raw iterator above. Its value type
	 * is an STL pair with the row and column indices, starting at 0, in the
	 * first and second positions, respectively.
	 */

	class RawIndexedIterator;
	class ConstRawIndexedIterator;

	RawIndexedIterator rawIndexedBegin ();
	RawIndexedIterator rawIndexedEnd ();
	ConstRawIndexedIterator rawIndexedBegin () const;
	ConstRawIndexedIterator rawIndexedEnd () const;

	/** Retrieve a row as a writeable reference
	 * @param i Row index
	 */
	Row &getRow (size_t i)
		{ return _A[i]; }

	/** Construct the transpose of this matrix and place it in the
	 * DenseRowsMatrix given
	 */
	DenseRowsMatrix &transpose (DenseRowsMatrix &AT) const;

	//@}

    protected:

	friend class SparseMatrixWriteHelper<Element, Row>;
	friend class SparseMatrixReadWriteHelper<Element, Row>;

	Rep               _A;
	size_t            _m;
	size_t            _n;
};

template <class Row>
struct MatrixTraits< DenseRowsMatrix<Row> >
{ 
	typedef DenseRowsMatrix<Row> MatrixType;
	typedef typename MatrixCategories::RowMatrixTag<MatrixTraits<MatrixType> > MatrixCategory; 
};

} // namespace LinBox

#endif // __LINBOX_dense_rows_matrix_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
