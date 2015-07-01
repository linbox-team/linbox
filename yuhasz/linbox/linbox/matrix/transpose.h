/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/matrix/transpose.h
 * Copyright (C) 2002 Bradford Hovinen,
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 *
 * Evolved from dense-base.h by Bradford Hovinen
 *
 * --------------------------------------------------------
 * 2003-01-11  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Move from blackbox/transpose-matrix.h to matrix/transpose.h
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __MATRIX_TRANSPOSE_H
#define __MATRIX_TRANSPOSE_H

#include <iostream>
#include <vector>
#include <fstream>

#include "linbox/vector/subiterator.h"
#include "linbox/vector/subvector.h"
#include "linbox/vector/stream.h"
#include "linbox/matrix/matrix-domain.h"

namespace LinBox
{

/** Matrix transpose
 * 
 * This class takes a matrix meeting the @ref{DenseMatrixBase} archetype and
 * switches the row and column iterators, giving the transpose of the original
 * matrix. It is generic with respect to the matrix given.
 * 
 * If the matrix given has limited iterators, then its transpose will have
 * limited iterators as well. In particular, if the matrix given has only row
 * iterators, then the transpose object will have only column iterators, and
 * vice versa.
 * 
 * This class differs from @ref{Transpose} in that it constructs a full matrix
 * representation, with row and/or column iterators. It does not include any
 * logic for matrix-vector products, and does not meet the
 * @ref{BlackboxArchetype} interface. Nor does it make such assumptions about
 * the matrix given.
 *
 * This class gives a constant matrix as output. It provides no iterators for
 * modification of the data in the matrix.
 * 
 * The input/output functionality of this class passes requests directly through
 * to the underlying matrix. In particular, the output will be the transpose of
 * the matrix expected and the input will expect the transpose of the matrix
 * given. Thus, it is not recommended to use TransposeMatrix for reading and
 * writing matrices, except for testing purposes.
 */
  
template <class Matrix, class Trait = typename MatrixTraits<Matrix>::MatrixCategory>
class TransposeMatrix
{
    public:

	typedef typename Matrix::Element Element;

	typedef typename Matrix::ConstColIterator ConstRowIterator;
	typedef typename Matrix::ConstRowIterator ConstColIterator;
	typedef typename Matrix::ConstRawIterator ConstRawIterator;
	typedef typename Matrix::ConstRawIndexedIterator ConstRawIndexedIterator;

	typedef typename Matrix::Row Column;
	typedef typename Matrix::Row Col;
	typedef typename Matrix::Column Row;

	/** Constructor.
	 * @param  A  Underlying matrix of which to construct the transpose
	 */
	TransposeMatrix (const Matrix &A)
		: _A (A)
	{}

	/** Copy constructor
	 */
	TransposeMatrix (const TransposeMatrix &M)
		: _A (M._A)
	{}

	/** Get the number of rows in the matrix
	 * @return Number of rows in matrix
	 */
	inline size_t rowdim () const
		{ return _A.coldim (); }

	/** Get the number of columns in the matrix
	 * @return Number of columns in matrix
	 */
	inline size_t coldim () const
		{ return _A.rowdim (); }

	/** @name Matrix I/O
	 */

	//@{

	/** Write a matrix to an output stream
	 * @param stream Stream to which to write the matrix
	 * @return Reference to stream
	 */
	template <class Field>
	inline std::ostream &write (std::ostream &stream, const Field &F) const
		{ return _A.write (stream, F); }

	//@} Matrix I/O

	/** @name Access to matrix elements
	 */

	//@{

	/** Get a read-only reference to the entry in the (i, j) position.
	 * @param i Row index
	 * @param j Column index
	 * @return Const reference to matrix entry
	 */
	inline const Element &getEntry (size_t i, size_t j) const
		{ return _A.getEntry (j, i); }

	/** Copy the (i, j) entry into x, and return a reference to x.
	 * This form is more in the Linbox style and is provided for interface
	 * compatibility with other parts of the library
	 * @param x Element in which to store result
	 * @param i Row index
	 * @param j Column index
	 * @return Reference to x
	 */
	inline Element &getEntry (Element &x, size_t i, size_t j) const
		{ return _A.getEntry (x, j, i); }

	/** @name Column of rows iterator
	 * The column of rows iterator traverses the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in dense format
	 */

	inline ConstRowIterator rowBegin () const { return _A.colBegin (); }
	inline ConstRowIterator rowEnd () const { return _A.colEnd (); }

	/** @name Row of columns iterator
	 * The row of columns iterator traverses the columns of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a column vector in dense format
	 */

	inline ConstColIterator colBegin () const { return _A.rowBegin (); }
	inline ConstColIterator colEnd () const { return _A.rowEnd (); }

	/** @name Raw iterator
	 *
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */

	inline ConstRawIterator rawBegin () const { return A.rawBegin (); }
	inline ConstRawIterator rawEnd () const { return A.rawEnd (); }

	/** @name Raw Indexed iterator
	 * Like the raw iterator, the indexed iterator is a method for 
	 * accessing all entries in the matrix in some unspecified order. 
	 * At each position of the the indexed iterator, it also provides 
	 * the row and column indices of the currently referenced entry.
	 * This is provided through it's rowIndex() and colIndex() functions.
	 */

	inline ConstRawIndexedIterator rawIndexedBegin() const { return A.rawIndexedBegin (); }
        inline ConstRawIndexedIterator rawIndexedEnd() const { return A.rawIndexedEnd (); }

	//@}

    protected:

	const Matrix &_A;
};

// Specialization for matrices that have both row and column iterators

template <class Matrix, class Trait>
class TransposeMatrix<Matrix, MatrixCategories::RowColMatrixTag<Trait> >
{
    public:

	typedef typename Matrix::Element Element;

	typedef typename Matrix::ConstColIterator ConstRowIterator;
	typedef typename Matrix::ConstRowIterator ConstColIterator;
	typedef typename Matrix::ConstRawIterator ConstRawIterator;
	typedef typename Matrix::ConstRawIndexedIterator ConstRawIndexedIterator;

	typedef typename Matrix::Row Column;
	typedef typename Matrix::Row Col;
	typedef typename Matrix::Column Row;

	TransposeMatrix (const Matrix &A) : _A (A) {}
	TransposeMatrix (const TransposeMatrix &M) : _A (M._A) {}

	inline size_t rowdim () const { return _A.coldim (); }
	inline size_t coldim () const { return _A.rowdim (); }

	template <class Field>
	inline std::ostream &write (std::ostream &stream, const Field &F) const
		{ return _A.write (stream, F); }

	inline const Element &getEntry (size_t i, size_t j) const { return _A.getEntry (j, i); }
	inline Element &getEntry (Element &x, size_t i, size_t j) const { return _A.getEntry (x, j, i); }

	inline ConstRowIterator rowBegin () const { return _A.colBegin (); }
	inline ConstRowIterator rowEnd () const { return _A.colEnd (); }

	inline ConstColIterator colBegin () const { return _A.rowBegin (); }
	inline ConstColIterator colEnd () const { return _A.rowEnd (); }

	inline ConstRawIterator rawBegin () const { return A.rawBegin (); }
	inline ConstRawIterator rawEnd () const { return A.rawEnd (); }

	inline ConstRawIndexedIterator rawIndexedBegin() const { return A.rawIndexedBegin (); }
        inline ConstRawIndexedIterator rawIndexedEnd() const { return A.rawIndexedEnd (); }

    protected:

	const Matrix &_A;
};

// Specialization for matrices that have only row iterators

template <class Matrix, class Trait>
class TransposeMatrix<Matrix, MatrixCategories::RowMatrixTag<Trait> >
{
    public:

	typedef typename Matrix::Element Element;

	typedef typename Matrix::ConstRowIterator ConstColIterator;
	typedef typename Matrix::ConstRawIterator ConstRawIterator;
	typedef typename Matrix::ConstRawIndexedIterator ConstRawIndexedIterator;

	typedef typename Matrix::Row Column;
	typedef typename Matrix::Row Col;

	TransposeMatrix (const Matrix &A) : _A (A) {}
	TransposeMatrix (const TransposeMatrix &M) : _A (M._A) {}

	inline size_t rowdim () const { return _A.coldim (); }
	inline size_t coldim () const { return _A.rowdim (); }

	template <class Field>
	inline std::ostream &write (std::ostream &stream, const Field &F) const
		{ return _A.write (stream, F); }

	inline const Element &getEntry (size_t i, size_t j) const { return _A.getEntry (j, i); }
	inline Element &getEntry (Element &x, size_t i, size_t j) const { return _A.getEntry (x, j, i); }

	inline ConstColIterator colBegin () const { return _A.rowBegin (); }
	inline ConstColIterator colEnd () const { return _A.rowEnd (); }

	inline ConstRawIterator rawBegin () const { return A.rawBegin (); }
	inline ConstRawIterator rawEnd () const { return A.rawEnd (); }

	inline ConstRawIndexedIterator rawIndexedBegin() const { return A.rawIndexedBegin (); }
        inline ConstRawIndexedIterator rawIndexedEnd() const { return A.rawIndexedEnd (); }

    protected:

	const Matrix &_A;
};

// Specialization for matrices that have only column iterators

template <class Matrix, class Trait>
class TransposeMatrix<Matrix, MatrixCategories::ColMatrixTag<Trait> >
{
    public:

	typedef typename Matrix::Element Element;

	typedef typename Matrix::ConstColIterator ConstRowIterator;
	typedef typename Matrix::ConstRawIterator ConstRawIterator;
	typedef typename Matrix::ConstRawIndexedIterator ConstRawIndexedIterator;

	typedef typename Matrix::Column Row;

	TransposeMatrix (const Matrix &A) : _A (A) {}
	TransposeMatrix (const TransposeMatrix &M) : _A (M._A) {}

	inline size_t rowdim () const { return _A.coldim (); }
	inline size_t coldim () const { return _A.rowdim (); }

	template <class Field>
	inline std::ostream &write (std::ostream &stream, const Field &F) const
		{ return _A.write (stream, F); }

	inline const Element &getEntry (size_t i, size_t j) const { return _A.getEntry (j, i); }
	inline Element &getEntry (Element &x, size_t i, size_t j) const { return _A.getEntry (x, j, i); }

	inline ConstRowIterator rowBegin () const { return _A.colBegin (); }
	inline ConstRowIterator rowEnd () const { return _A.colEnd (); }

	inline ConstRawIterator rawBegin () const { return A.rawBegin (); }
	inline ConstRawIterator rawEnd () const { return A.rawEnd (); }

	inline ConstRawIndexedIterator rawIndexedBegin() const { return A.rawIndexedBegin (); }
        inline ConstRawIndexedIterator rawIndexedEnd() const { return A.rawIndexedEnd (); }

    protected:

	const Matrix &_A;
};

template <class Matrix, class Trait>
struct MatrixTraits< TransposeMatrix<Matrix, MatrixCategories::RowColMatrixTag<Trait> > >
{ 
	typedef TransposeMatrix<Matrix, MatrixCategories::RowColMatrixTag<Trait> > MatrixType;
	typedef typename MatrixCategories::RowColMatrixTag<MatrixTraits<MatrixType> > MatrixCategory; 
};

template <class Matrix, class Trait>
struct MatrixTraits< TransposeMatrix<Matrix, MatrixCategories::RowMatrixTag<Trait> > >
{ 
	typedef TransposeMatrix<Matrix, MatrixCategories::RowMatrixTag<Trait> > MatrixType;
	typedef typename MatrixCategories::ColMatrixTag<MatrixTraits<MatrixType> > MatrixCategory; 
};

template <class Matrix, class Trait>
struct MatrixTraits< TransposeMatrix<Matrix, MatrixCategories::ColMatrixTag<Trait> > >
{ 
	typedef TransposeMatrix<Matrix, MatrixCategories::ColMatrixTag<Trait> > MatrixType;
	typedef typename MatrixCategories::RowMatrixTag<MatrixTraits<MatrixType> > MatrixCategory; 
};

} // namespace LinBox

#endif // __MATRIX_TRANSPOSE_INL
