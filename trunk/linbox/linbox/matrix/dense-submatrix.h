/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/dense-submatrix.h
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
 * -----------------------------------------------------------
 * 2002-10-27  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 * 
 * Rename from densesubmatrix.h
 *
 * Constructor modifications: changed the interface to match Submatrix
 * -----------------------------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __DENSE_SUBMATRIX_H
#define __DENSE_SUBMATRIX_H

#include "linbox-config.h"

#include "linbox/util/debug.h"
#include "linbox/blackbox/dense.h"
#include "linbox/blackbox/archetype.h"

namespace LinBox
{

/** Submatrix of a dense matrix
 *
 * This matrix type conforms to the same interface as @ref{DenseMatrix}, except
 * that you cannot resize it. It represents a submatrix of a dense matrix. Upon
 * construction, one can freely manipulate the entries in the DenseSubmatrix,
 * and the corresponding entries in the underlying DenseMatrix will be modified.
 */
template<class Field, class Vector = typename LinBox::Vector<Field>::Dense>
class DenseSubmatrix : public BlackboxArchetype<Vector>
{
    public:
 
	typedef typename Field::Element Element;
	typedef typename Vector::iterator pointer;

	class RawIterator;
	class ConstRawIterator;

	typedef typename DenseMatrix<Field>::RowIterator            RowIterator;
	typedef typename DenseMatrix<Field>::ConstRowIterator       ConstRowIterator;
	typedef typename DenseMatrix<Field>::Row                    Row;
	typedef typename DenseMatrix<Field>::ConstRow               ConstRow;
	typedef typename DenseMatrix<Field>::ColOfRowsIterator      ColOfRowsIterator;
	typedef typename DenseMatrix<Field>::ConstColOfRowsIterator ConstColOfRowsIterator;
	typedef typename DenseMatrix<Field>::ColIterator            ColIterator;
	typedef typename DenseMatrix<Field>::ConstColIterator       ConstColIterator;
	typedef typename DenseMatrix<Field>::Col                    Col;
	typedef typename DenseMatrix<Field>::ConstCol               ConstCol;
	typedef typename DenseMatrix<Field>::RowOfColsIterator      RowOfColsIterator;
	typedef typename DenseMatrix<Field>::ConstRowOfColsIterator ConstRowOfColsIterator;

	/** Empty constructor
	 */
	DenseSubmatrix () {}

	/** Constructor from an existing @ref{DenseMatrix} and dimensions
	 * @param M Pointer to @ref{DenseMatrix} of which to construct submatrix
	 * @param row Starting row
	 * @param col Starting column
	 * @param rowdim Row dimension
	 * @param coldim Column dimension
	 */
	DenseSubmatrix (DenseMatrix<Field, Vector> *M,
			size_t row,
			size_t col,
			size_t rowdim,
			size_t coldim);

	/** Constructor from an existing submatrix and dimensions
	 * @param SM Constant reference to DenseSubmatrix from which to
	 *           construct submatrix
	 * @param row Starting row
	 * @param col Starting column
	 * @param rowdim Row dimension
	 * @param coldim Column dimension
	 */
	DenseSubmatrix (const DenseSubmatrix<Field, Vector> &SM,
			size_t row,
			size_t col,
			size_t rowdim,
			size_t coldim);

	/** Copy constructor
	 * @param _M Submatrix to copy
	 */
	DenseSubmatrix (const DenseSubmatrix<Field, Vector> &SM);

	/** Assignment operator
	 * Assign the given submatrix to this one
	 * @param _M Submatrix to assign
	 * @return Reference to this submatrix
	 */
	DenseSubmatrix &operator = (const DenseSubmatrix<Field, Vector> &SM);

	/** Construct a clone of the submatrix
	 * @return Pointer to clone of the submatrix
	 */
	BlackboxArchetype<Vector>* clone () const
		{ return new DenseSubmatrix<Field, Vector> (*this); }
       
	/** Get the number of rows in the matrix
	 * @return Number of rows in matrix
	 */
	size_t rowdim () const
		{ return _end_row - _beg_row; }

	/** Get the number of columns in the matrix
	 * @return Number of columns in matrix
	 */
	size_t coldim () const
		{ return _end_col - _beg_col; }

	/** Retrieve the field over which this matrix is defined
	 * @return Reference to the underlying field
	 */
	const Field &field () const
		{ return _M->field ();}
    
	/** @name Input and output
	 */

	//@{

	/** Read the matrix from an input stream
	 * @param file Input stream from which to read
	 */
	void read (std::istream &file);
    
	/** Write the matrix to an output stream
	 * @param os Output stream to which to write
	 */
	std::ostream &write (std::ostream &os = std::cout) const;

	//@}

	/** @name Access to matrix elements
	 */

	//@{
    
	/** Set the entry at (i, j)
	 * @param i Row number, 0...rowdim () - 1
	 * @param j Column number 0...coldim () - 1
	 * @param a_ij Element to set
	 */
	void setEntry (size_t i, size_t j, const Element &a_ij)
		{ _M->setEntry (_beg_row + i, _beg_col + j, a_ij); }

	/** Get a writeable reference to an entry in the matrix
	 * @param i Row index of entry
	 * @param j Column index of entry
	 * @return Reference to matrix entry
	 */
	Element &refEntry (size_t i, size_t j)
		{ return _M->refEntry (i + _beg_row, j + _beg_col); } 

	/** Get a read-only individual entry from the matrix
	 * @param i Row index
	 * @param j Column index
	 * @return Const reference to matrix entry
	 */
	const Element &getEntry (size_t i, size_t j) const
		{ return _M->getEntry (i + _beg_row, j + _beg_col); } 

	/** Get an entry and store it in the given value
	 * This form is more in the Linbox style and is provided for interface
	 * compatibility with other parts of the library
	 * @param x Element in which to store result
	 * @param i Row index
	 * @param j Column index
	 * @return Reference to x
	 */
	Element &getEntry (Element &x, size_t i, size_t j)
		{ return _M->getEntry (a, i + _beg_row, j + _beg_col); } 

	/** @name Columns of rows iterator
	 * The columns of row iterator gives each of the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in dense format
	 */

	ColOfRowsIterator colOfRowsBegin ();
	ColOfRowsIterator colOfRowsEnd ();
	ConstColOfRowsIterator colOfRowsBegin () const;
	ConstColOfRowsIterator colOfRowsEnd () const;
 
	/** @name Row of columns iterator
	 * The row of columns iterator gives each of the columns of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a column vector in dense format
	 */

	RowOfColsIterator rowOfColsBegin ();
	RowOfColsIterator rowOfColsEnd ();
	ConstRowOfColsIterator rowOfColsBegin () const;
	ConstRowOfColsIterator rowOfColsEnd () const;

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

	/** @name Index iterator
	 * The index iterator gives the row, column indices of all matrix
	 * elements in the same order as the raw iterator above. Its value type
	 * is an STL pair with the row and column indices, starting at 0, in the
	 * first and second positions, respectively.
	 */

        class RawIndexIterator;
        typedef const RawIndexIterator ConstRawIndexIterator;

        RawIndexIterator rawIndexBegin();
        RawIndexIterator rawIndexEnd();   
	ConstRawIndexIterator rawIndexBegin() const;
        ConstRawIndexIterator rawIndexEnd() const;   

	/** Retrieve a reference to a row
	 * @param i Row index
	 */
	Row operator[] (int i);
	ConstRow operator[] (int i) const;

	//@}

	/** @name Black box interface
	 */

	//@{

	/** Generic matrix-vector apply
	 * y = A * x
	 * This version of apply allows use of arbitrary input and output vector
	 * types
	 * @param y Output vector
	 * @param x Input vector
	 * @return Reference to output vector
	 */
	template<class Vect1, class Vect2>
	Vect1 &apply (Vect1 &y, const Vect2 &x) const;

	/** Generic in-place apply
	 * y = A * y
	 * This version of in-place apply allows use of an arbitrary vector
	 * type. Because it performs allocation and copying, it is not
	 * recommended for general use.
	 * @param y Input vector
	 * @return Reference to output vector
	 */
	template<class Vect1>
	Vect1 &applyIn (Vect1 &y) const
	{
		std::vector<Element> x (y.begin (), y.end ());
		apply (y, x);
		return y;
	}

	/** Matrix-vector apply
	 * y = A * x
	 * This implements the @ref{BlackboxArchetype} apply requirement
	 * @param y Output vector
	 * @param x Input vector
	 * @return Reference to output vector
	 */
	Vector &apply (Vector &y, const Vector &x) const
		{ return apply<Vector,Vector> (y, x); }

	/** Iterator form of apply
	 * This form of apply takes iterators specifying the beginning and end
	 * of the vector to which to apply the matrix, and the beginning of the
	 * vector at which to store the result of application. It is generic
	 * with respect to iterator type, allowing different iterators to be
	 * used for the input and output vectors.
	 * @param out Beginning of output vector
	 * @param inbegin Beginning of input vector
	 * @param outbegin End of input vector
	 * @return Reference to beginning of output vector
	 */
	template<class Iterator1, class Iterator2>
	Iterator1 &apply (Iterator1 in, const Iterator2 &outbegin, const Iterator2 &outend) const;

	/** Generic matrix-vector transpose apply
	 * y = A^T * x
	 * This version of applyTranspose allows use of arbitrary input and
	 * output vector types
	 * @param y Output vector
	 * @param x Input vector
	 * @return Reference to output vector
	 */
	template<class Vect1, class Vect2>
	Vect1 &applyTranspose (Vect1 &y, const Vect2 &x) const;

	/** Generic in-place transpose apply
	 * y = A^T * y
	 * This version of in-place transpose apply allows use of an arbitrary
	 * vector type. Because it performs allocation and copying, it is not
	 * recommended for general use.
	 * @param y Input vector
	 * @return Reference to output vector
	 */
	template<class Vect>
	Vect &applyTransposeIn (Vect &y) const
	{
		std::vector<Element> x (y.begin (), y.end ());
		applyTranspose (y, x);
		return y;
	}

	/** Matrix-vector transpose apply
	 * y = A^T * x
	 * This implements the @ref{BlackboxArchetype} applyTranspose
	 * requirement
	 * @param y Output vector
	 * @param x Input vector
	 * @return Reference to output vector
	 */
	Vector &applyTranspose (Vector &y, const Vector &x) const
		{ return applyTranspose<Vector,Vector> (y, x); }

	/** Iterator form of transpose apply
	 *
	 * This form of transpose apply takes iterators specifying the beginning
	 * and end of the vector to which to apply the matrix, and the beginning
	 * of the vector at which to store the result of application. It is
	 * generic with respect to iterator type, allowing different iterators
	 * to be used for the input and output vectors.
	 *
	 * @param out Beginning of output vector
	 * @param inbegin Beginning of input vector
	 * @param outbegin End of input vector
	 * @return Reference to beginning of output vector
	 */
	template<class Iterator1, class Iterator2>
	Iterator1 &applyTranspose (Iterator1 in, const Iterator2 &outbegin, const Iterator2 &outend) const;

	//@}

    protected:
	DenseMatrix<Field> *_M;
	size_t _beg_row;
	size_t _end_row;
	size_t _beg_col;
	size_t _end_col;
};

} // namespace LinBox

#include "linbox/blackbox/dense-submatrix.inl"

#endif // __DENSE_SUBMATRIX_H

