/* linbox/matrix/archetype.h
 * Copyright (C) 2001 B. David Saunders,
 *               2001-2002 Bradford Hovinen,
 *               2002 Zhendong Wan
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>,
 *            Zhendong Wan <wan@mail.eecis.udel.edu>
 *
 * Borrowed from dense-base.h by Bradford Hovinen
 * evolved from dense-matrix.h by -bds, Zhendong Wan
 *
 * This holds the "directly represented" matrix archetype. It is provided here
 * only for reference; it does not provide any useful functionality. See the
 * other headers in this directory for useful classes.
 *
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __LINBOX_matrix_archetype_H
#define __LINBOX_matrix_archetype_H

#include <iostream>
#include <vector>
#include <fstream>

#include "linbox/blackbox/archetype.h"

namespace LinBox
{

/** @brief Directly-represented matrix archetype
 *
 * This archetype gives the common interface for matrices that have direct
 * representations. The matrices are required to provide iterators to access and
 * manipulate their entries, but not any matrix-matrix or matrix-vector
 * arithmetic. That is, they are pure containers. As such, they are only
 * parameterized on the element type, not on the field type.
 */
  
template <class _Element>
class MatrixArchetype
{
    public:

	typedef _Element Element;

	/** Empty Constructor.
	 */
	MatrixArchetype ();

	/** Constructor with size
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	MatrixArchetype (size_t m, size_t n);

	/** Copy constructor
	 */
	MatrixArchetype (const MatrixArchetype &M);

	/** Operator =
	 */
	MatrixArchetype& operator= (const MatrixArchetype& M);

	/** Get the number of rows in the matrix
	 * @return Number of rows in matrix
	 */
	size_t rowdim () const;

	/** Get the number of columns in the matrix
	 * @return Number of columns in matrix
	 */
	size_t coldim () const;

	/** \brief Resize the matrix to the given dimensions
	 *
	 * The state of the matrix's entries after a call to this method is
	 * undefined.
	 *
	 * This interface is optional; a matrix can omit it if it makes no sense
	 * in the context.
	 *
	 * @param m Number of rows
	 * @param n Number of columns
	 */
	void resize (size_t m, size_t n);

	/** @name Input and output
	 */

	//@{

	/** Read the matrix from an input stream
	 * @param file Input stream from which to read
	 * @param F Field over which to read
	 */
	template <class Field>
	std::istream &read (std::istream &file, const Field &F);

	/** Write the matrix to an output stream
	 * @param os Output stream to which to write
	 * @param F Field over which to write
	 */
	template <class Field>
	std::ostream &write (std::ostream &os, const Field &F) const;

	//@}

	/** @name Access to matrix elements
	 */

	//@{

	/** Set the entry at the (i, j) position to a_ij.
	 * @param i Row number, 0...rowdim () - 1
	 * @param j Column number 0...coldim () - 1
	 * @param a_ij Element to set
	 */
	void setEntry (size_t i, size_t j, const Element &a_ij);

	/** Get a writeable reference to the entry in the (i, j) position.
	 * @param i Row index of entry
	 * @param j Column index of entry
	 * @return Reference to matrix entry
	 */
	Element &refEntry (size_t i, size_t j);

	/** Get a read-only reference to the entry in the (i, j) position.
	 * @param i Row index
	 * @param j Column index
	 * @return Const reference to matrix entry
	 */
	const Element &getEntry (size_t i, size_t j) const;

	/** Copy the (i, j) entry into x, and return a reference to x.
	 * This form is more in the Linbox style and is provided for interface
	 * compatibility with other parts of the library
	 * @param x Element in which to store result
	 * @param i Row index
	 * @param j Column index
	 * @return Reference to x
	 */
	Element &getEntry (Element &x, size_t i, size_t j) const;

	/* N.B. A matrix type may omit either one, but not both, of the
	 * following two iterator types. If one type is omitted, then certain
	 * restrictions on matrix-matrix arithmetic apply; see
	 * @ref{MatrixDomain}
	 */

	/** @name Column of rows iterator
	 * The column of rows iterator traverses the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in dense format
	 */

	class Row;
	class ConstRow;  
	class RowIterator;    
	class ConstRowIterator;

	RowIterator rowBegin ();  
	RowIterator rowEnd ();
	ConstRowIterator rowBegin () const;        
	ConstRowIterator rowEnd () const;

	/** @name Row of columns iterator
	 * The row of columns iterator traverses the columns of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a column vector in dense format
	 */

	class Col;
	class ConstCol;
	class ColIterator;
	class ConstColIterator;
    
	typedef Col Column;
	typedef ConstCol ConstColumn;

	ColIterator colBegin ();
	ColIterator colEnd ();
	ConstColIterator colBegin () const;    
	ConstColIterator colEnd () const;

	/** @name Raw iterator
	 *
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */

	class RawIterator;
	class ConstRawIterator;
    
	RawIterator rawBegin ();		  
	RawIterator rawEnd ();
	ConstRawIterator rawBegin () const;
	ConstRawIterator rawEnd () const;

	/** @name Raw Indexed iterator
	 * Like the raw iterator, the indexed iterator is a method for 
	 * accessing all entries in the matrix in some unspecified order. 
	 * At each position of the the indexed iterator, it also provides 
	 * the row and column indices of the currently referenced entry.
	 * This is provided through it's rowIndex() and colIndex() functions.
	 */

        class RawIndexedIterator;
        typedef const RawIndexedIterator ConstRawIndexedIterator;

        RawIndexedIterator rawIndexedBegin();
        RawIndexedIterator rawIndexedEnd();   
	ConstRawIndexedIterator rawIndexedBegin() const;
        ConstRawIndexedIterator rawIndexedEnd() const;   
    
	/** Retrieve a reference to a row.
	 * Since rows may also be indexed, this allows A[i][j] notation
	 * to be used.
	 *
	 * This may be omitted by an implementation if no Row type is available
	 *
	 * @param i Row index
	 */
	Row operator[] (size_t i);
	ConstRow operator[] (size_t i) const;

	//@}

	/** @name Computing matrix information
	 */

	//@{

	/** Compute the column density, i.e. the number of entries per column
	 */
	template <class Vector>
	Vector &columnDensity (Vector &v) const;

	/** Compute the transpose
	 */
	MatrixArchetype &transpose (MatrixArchetype &M) const;

	//@}

    protected:

	std::vector<Element>  _rep;
	size_t                _rows, _cols;
};

template <class Element>
struct MatrixTraits< MatrixArchetype<Element> >
{ 
	typedef MatrixArchetype<Element> MatrixType;
	typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
};

} // namespace LinBox

#endif // __LINBOX_matrix_archetype_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
