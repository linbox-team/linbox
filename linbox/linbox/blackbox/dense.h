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
  
template <class Field>
class DenseMatrix : public BlackboxArchetype<std::vector<typename Field::Element> >
{
    public:
	typedef typename Field::Element                 Element;
	typedef typename std::vector<Element>           Vector;
	typedef typename std::vector<Element>::iterator pointer;
    
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
	DenseMatrix (const Field &F, size_t m, size_t n, RandIter& iter)
		: _F (F), _VD (F), _rep (m * n), _rows (m), _cols (n)
	{
		for (typename Vector::iterator p = _rep.begin (); p != _rep.end (); ++p)
			iter.random (*p);
	}
    
	/** Constructor.
	 * @param  F the field of entries; passed so that arithmetic may be done on elements. 
	 * @param  stream A vector stream to use as a source of vectors for this matrix
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
    
	/// Blackbox interface
    
	BlackboxArchetype<Vector> *clone () const 
		{ return new DenseMatrix<Field> (*this); }
    
	template<class Vect1, class Vect2>
	Vect1& apply (Vect1& y, const Vect2& x) const;
    
	Vector& apply (Vector &y, const Vector &x) const
		{ return apply<Vector,Vector> (y, x); }
    
	template<class Vect1>
	Vect1& applyIn (Vect1& y) const
	{
		std::vector<Element> x (y.begin (),y.end ());
		apply (y,x);
		return y;
	}
      
	Vector& applyIn (Vector& y) const
	{ return applyIn<Vector> (y);}

	template<class Iterator1, class Iterator2 >
	Iterator1& apply (Iterator1 out, 
			  const Iterator2& inbegin, 
			  const Iterator2& inend) const;
   
	template<class Vect1, class Vect2>
	Vect1& applyTranspose (Vect1& y, const Vect2& x) const;
    
	Vector& applyTranspose (Vector& y, const Vector& x) const
		{ return applyTranspose<Vector,Vector> (y, x); }

	template<class Vect>
	Vect& applyTransposeIn (Vect& y) const
	{
		std::vector<Element> x (y.begin (), y.end ());
		applyTranspose (y, x);
		return y;
	}
  
	Vector& applyTransposeIn (Vector& y) const
		{ return applyTransposeIn<Vector> (y); }
    
	template<class Iterator1, class Iterator2>
	Iterator1& applyTranspose (Iterator1 out, 
				   const Iterator2& inbegin, 
				   const Iterator2& inend) const;

	size_t rowdim (void) const;
	size_t coldim (void) const;
    
	/// entry access raw view.  Size m*x vector in C (row major) order.
	typedef typename Vector::iterator RawIterator;
	typedef typename Vector::const_iterator ConstRawIterator;
    
	RawIterator rawBegin ();		  
	RawIterator rawEnd ();
    
	ConstRawIterator rawBegin () const;
	ConstRawIterator rawEnd () const;

	/** @name Index iterator
         * The index iterator gives the row, column indices of all matrix
         * elements in the same order as the raw iterator above.
         */
        class RawIndexIterator;
        typedef const RawIndexIterator ConstRawIndexIterator;

        RawIndexIterator rawIndexBegin();
        RawIndexIterator rawIndexEnd();   
	ConstRawIndexIterator rawIndexBegin() const;
        ConstRawIndexIterator rawIndexEnd() const;   
  
	// col sequence of rows view
	typedef typename Vector::iterator RowIterator;
	typedef typename Vector::const_iterator ConstRowIterator;
    
	class ColOfRowsIterator;
    
	typedef Subvector<DenseMatrix<Field>::RowIterator> Row;  
	typedef Subvector<DenseMatrix<Field>::ConstRowIterator> ConstRow;
     
	class ColOfRowsIterator;    
	class ConstColOfRowsIterator;

	ColOfRowsIterator colOfRowsBegin ();  
	ColOfRowsIterator colOfRowsEnd ();
 
	ConstColOfRowsIterator colOfRowsBegin () const;        
	ConstColOfRowsIterator colOfRowsEnd () const;
    

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
    
	/** Set the entry at (i, j)
	 * @param i Row number, 0...rowdim () - 1
	 * @param j Column number 0...coldim () - 1
	 * @param a_ij Element to set
	 */
	void setEntry (size_t i, size_t j, const Element& a_ij) ;
    
	Element& getEntry (Element& a_ij, size_t i, size_t j) const;
    
	/** Read the matrix from an input stream
	 * @param file Input stream from which to read
	 */
	void read (std::istream &file);
    
	/** Write the matrix to an output stream
	 * @param os Output stream to which to write
	 */
	std::ostream &write (std::ostream &os = std::cout) const;
 
	const Field& field () const
		{ return _F;}
    
	Row operator[] (size_t i);
	ConstRow operator[] (size_t i) const;
   
    protected:
    
	const Field          &_F;
	VectorDomain<Field>   _VD;
	std::vector<Element>  _rep;
	size_t                _rows, _cols;
};

}

#include "dense.inl"

#endif
