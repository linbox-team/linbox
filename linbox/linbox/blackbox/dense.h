/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/dense.h
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
 * 2002-08-09  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Renamed file from dense-matrix1.h to dense.h
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __DENSE_H
#define __DENSE_H

#include <iostream>
#include <vector>
#include <fstream>

#include "linbox/blackbox/archetype.h"
#include "linbox/vector/subiterator.h"
#include "linbox/vector/subvector.h"
#include "linbox/vector/stream.h"
#include "linbox/field/vector-domain.h"

namespace LinBox
{

/** Blackbox dense matrix template. This is a class of dense matrices
 * templatized by the {@link Fields field} in which the elements
 * reside. The matrix is stored as a one dimensional STL vector of
 * the elements, by rows. The interface provides for iteration
 * over rows and over columns.
 *
 * The class also conforms to the {@link Archetypes archetype} for
 * \Ref{Blackbox Matrices}.
 *
 * @param Field \Ref{LinBox} field
 */
  
template <class Field, class Vector = typename LinBox::Vector<Field>::Dense>
class DenseMatrix : public BlackboxArchetype<Vector>
{
    public:
	typedef typename Field::Element   Element;
	typedef typename Vector::iterator pointer;

	/** Constructor.
	 * @param  F the field of entries; passed so that arithmetic may be done on elements. 
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	DenseMatrix (const Field &F, size_t m, size_t n)
		: _F (F), _VD (F), _rep (m * n), _rows (m), _cols (n)
	{}
    
	/** Constructor.
	 * @param  F the field of entries; passed so that arithmetic may be done on elements. 
	 * @param  m  row dimension
	 * @param  n  column dimension
	 * @para iter, random iterator
	 */
	template<class RandIter>
	DenseMatrix (const Field &F, size_t m, size_t n, RandIter &iter)
		: _F (F), _VD (F), _rep (m * n), _rows (m), _cols (n)
	{
		for (typename Vector::iterator p = _rep.begin (); p != _rep.end (); ++p)
			iter.random (*p);
	}
    
	/** Constructor using a vector stream
	 * @param  F The field of entries; passed so that arithmetic may be done
	 *           on elements. 
	 * @param  stream A vector stream to use as a source of vectors for this
	 *                matrix
	 */
	template <class StreamVector>
	DenseMatrix (const Field &F, VectorStream<StreamVector> &stream)
		: _F (F), _VD (F), _rep (stream.dim () * stream.size ()), _rows (stream.size ()), _cols (stream.dim ())
	{
		StreamVector tmp;

		VectorWrapper::ensureDim (tmp, stream.dim ());

		for (ColOfRowsIterator p = colOfRowsBegin (); p != colOfRowsEnd (); ++p) {
			stream >> tmp;
			_VD.copy (*p, tmp);
		}
	}

	/** Copy constructor
	 */
	DenseMatrix (const DenseMatrix &M)
		: _F (M._F), _VD (M._F), _rep (M._rep),_rows (M._rows), _cols (M._cols)
	{}

	/** Construct a copy of the matrix and return a pointer to it
	 * @return Pointer to copy of the matrix
	 */
	BlackboxArchetype<Vector> *clone () const 
		{ return new DenseMatrix<Field, Vector> (*this); }

	/** Get the number of rows in the matrix
	 * @return Number of rows in matrix
	 */
	size_t rowdim () const
		{ return _rows; }

	/** Get the number of columns in the matrix
	 * @return Number of columns in matrix
	 */
	size_t coldim () const
		{ return _cols; }

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
		{ _rep[i*_cols+j] = a_ij; }

	/** Get a writeable reference to an entry in the matrix
	 * @param i Row index of entry
	 * @param j Column index of entry
	 * @return Reference to matrix entry
	 */
	Element &refEntry (size_t i, size_t j)
		{ return _rep[i * _cols + j]; }

	/** Get a read-only individual entry from the matrix
	 * @param i Row index
	 * @param j Column index
	 * @return Const reference to matrix entry
	 */
	const Element &getEntry (size_t i, size_t j) const
		{ return _rep[i * _cols + j]; }

	/** Get an entry and store it in the given value
	 * This form is more in the Linbox style and is provided for interface
	 * compatibility with other parts of the library
	 * @param x Element in which to store result
	 * @param i Row index
	 * @param j Column index
	 * @return Reference to x
	 */
	Element &getEntry (Element &x, size_t i, size_t j) const
		{ x = _rep[i * _cols + j]; return x; }

	/** @name Columns of rows iterator
	 * The columns of row iterator gives each of the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in dense format
	 */

	typedef typename Vector::iterator RowIterator;
	typedef typename Vector::const_iterator ConstRowIterator;

	typedef Subvector<RowIterator> Row;  
	typedef Subvector<ConstRowIterator> ConstRow;

	class ColOfRowsIterator;    
	class ConstColOfRowsIterator;

	ColOfRowsIterator colOfRowsBegin ();  
	ColOfRowsIterator colOfRowsEnd ();
	ConstColOfRowsIterator colOfRowsBegin () const;        
	ConstColOfRowsIterator colOfRowsEnd () const;

	/** @name Row of columns iterator
	 * The row of columns iterator gives each of the columns of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a column vector in dense format
	 */

	typedef Subiterator<typename Vector::iterator> ColIterator;
	typedef Subiterator<typename Vector::const_iterator> ConstColIterator;

	typedef Subvector<ColIterator> Col;
	typedef Subvector<ConstColIterator> ConstCol;
    
	class RowOfColsIterator;
	class ConstRowOfColsIterator;
    
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

	typedef typename Vector::iterator RawIterator;
	typedef typename Vector::const_iterator ConstRawIterator;
    
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
	Row operator[] (size_t i);
	ConstRow operator[] (size_t i) const;

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
		std::vector<Element> x (y.begin (),y.end ());
		apply (y,x);
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
		{ return apply<Vector, Vector> (y, x); }

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
	template<class Iterator1, class Iterator2 >
	Iterator1 &apply (Iterator1 out, 
			  const Iterator2 &inbegin, 
			  const Iterator2 &inend) const;

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
	Iterator1 &applyTranspose (Iterator1 out, 
				   const Iterator2 &inbegin, 
				   const Iterator2 &inend) const;

	//@}

	/** Retrieve the field over which this matrix is defined
	 * @return Reference to the underlying field
	 */
	const Field &field () const
		{ return _F;}

    protected:

	const Field          &_F;
	VectorDomain<Field>   _VD;
	std::vector<Element>  _rep;
	size_t                _rows, _cols;
};

}

#include "dense.inl"

#endif
