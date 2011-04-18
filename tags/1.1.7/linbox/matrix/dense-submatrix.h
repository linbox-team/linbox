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

#ifndef __LINBOX_dense_submatrix_H
#define __LINBOX_dense_submatrix_H

#include "linbox/linbox-config.h"

#include "linbox/util/debug.h"
#include "linbox/matrix/dense.h"
#include "linbox/matrix/matrix-domain.h"

namespace LinBox
{

/** Submatrix of a dense matrix
 *
 * This matrix type conforms to the same interface as @ref{DenseMatrixBase},
 * except that you cannot resize it. It represents a submatrix of a dense
 * matrix. Upon construction, one can freely manipulate the entries in the
 * DenseSubmatrix, and the corresponding entries in the underlying
 * DenseMatrixBase will be modified.

\ingroup matrix
 */
template<class _Element>
class DenseSubmatrix
{
    public:
 
		class RawIterator ;
		class ConstRawIterator;

		typedef _Element Element;
		typedef DenseSubmatrix<_Element> Self_t;
    

	/** \brief 
	 *
	 * The row iterator gives the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in dense format
	 */
	typedef typename DenseMatrixBase<Element>::RowIterator            RowIterator;
	typedef typename DenseMatrixBase<Element>::ConstRowIterator       ConstRowIterator;
	typedef typename DenseMatrixBase<Element>::Row                    Row;
	typedef typename DenseMatrixBase<Element>::ConstRow               ConstRow;

	/** \brief
	 *
	 * The columns iterator gives the columns of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a column vector in dense format
	 */
	typedef typename DenseMatrixBase<Element>::ColIterator            ColIterator;
	typedef typename DenseMatrixBase<Element>::ConstColIterator       ConstColIterator;
	typedef typename DenseMatrixBase<Element>::Col                    Col;
	typedef typename DenseMatrixBase<Element>::Column                 Column;
	typedef typename DenseMatrixBase<Element>::ConstCol               ConstCol;

	/** \brief
	 */
	DenseSubmatrix () :_M(NULL) {}

	/** Constructor from an existing @ref{DenseMatrixBase} and dimensions
	 * \param M Pointer to @ref{DenseMatrixBase} of which to construct submatrix
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

	/** Constructor from an existing @ref{DenseMatrixBase} 
	 * \param M Pointer to @ref{DenseMatrixBase} of which to construct submatrix
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

	/** Copy constructor
	 * @param _M Submatrix to copy
	 */
	DenseSubmatrix (const DenseSubmatrix<Element> &SM);

	/** Assignment operator
	 * Assign the given submatrix to this one
	 * @param _M Submatrix to assign
	 * @return Reference to this submatrix
	 */
	DenseSubmatrix &operator = (const DenseSubmatrix<Element> &SM);

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

	RawIterator rawBegin ();     
	RawIterator rawEnd ();
	ConstRawIterator rawBegin () const;       
	ConstRawIterator rawEnd () const;  


	template<typename _Tp1>
        struct rebind
        { 
            typedef DenseSubmatrix<typename _Tp1::Element> other;
        };




	/** Read the matrix from an input stream
	 * @param file Input stream from which to read
	 * @param field 
	 */
	template<class Field>
	std::istream& read (std::istream &file, const Field& field);
    
	/** Write the matrix to an output stream
	 * @param os Output stream to which to write
	 * @param field
	 */
	template<class Field>
	std::ostream& write (std::ostream &os, const Field& field, bool mapleFormat = false) const;
	
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
	Element &getEntry (Element &x, size_t i, size_t j) const
		{ return _M->getEntry (x, i + _beg_row, j + _beg_col); } 

	RowIterator rowBegin ();
	RowIterator rowEnd ();
	ConstRowIterator rowBegin () const;
	ConstRowIterator rowEnd () const;
 
	ColIterator colBegin ();
	ColIterator colEnd ();
	ConstColIterator colBegin () const;
	ConstColIterator colEnd () const;

	/** \brief
	 *
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */

	//class RawIterator;   
	//class ConstRawIterator;
   
	/** \brief
	 *
	 * Like the raw iterator, the indexed iterator is a method for
	 * accessing all entries in the matrix in some unspecified order.
	 * At each position of the the indexed iterator, it also provides
	 * the row and column indices of the currently referenced entry.
	 * This is provided through it's rowIndex() and colIndex() functions.
	 */

        class RawIndexedIterator;
        class ConstRawIndexedIterator;

        RawIndexedIterator rawIndexedBegin();
        RawIndexedIterator rawIndexedEnd();   
	ConstRawIndexedIterator rawIndexedBegin() const;
        ConstRawIndexedIterator rawIndexedEnd() const;   

	/** Retrieve a reference to a row
	 * @param i Row index
	 */
	//Row operator[] (int i);               not actually used, causes a compile error...
	//ConstRow operator[] (int i) const;

    protected:
	DenseMatrixBase<Element> *_M;
	size_t _beg_row;
	size_t _end_row;
	size_t _beg_col;
	size_t _end_col;
};

template <class Element>
struct MatrixTraits< DenseSubmatrix<Element> >
{ 
	typedef DenseSubmatrix<Element> MatrixType;
	typedef typename MatrixCategories::RowColMatrixTag MatrixCategory; 
};

} // namespace LinBox

#include "linbox/matrix/dense-submatrix.inl"

#endif // __LINBOX_dense_submatrix_H


/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
