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
 * 2002-10-27  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Split out container/iterator functionality into DenseMatrixBase
 * --------------------------------------------------------
 * 2002-08-09  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Renamed file from dense-matrix1.h to dense.h
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __BLACKBOX_DENSE_H
#define __BLACKBOX_DENSE_H

#include <iostream>
#include <vector>
#include <fstream>

#include "linbox/vector/subiterator.h"
#include "linbox/vector/subvector.h"
#include "linbox/vector/stream.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/matrix/dense.h"
#include <linbox/matrix/matrix-domain.h>

namespace LinBox
{

/** @memo Blackbox interface to dense matrix representation. 
 * @doc This is a class of dense matrices
 * templatized by the {@link Fields field} in which the elements
 * reside. The matrix is stored as a one dimensional STL vector of
 * the elements, in row major order. The interface provides for iteration
 * over rows and over columns.
 *
 * The class also conforms to the {@link Archetypes archetype} for
 * \Ref{Blackbox Matrices}.
 *
 * Currently, only dense vectors are supported when doing matrix-vector
 * applies.
 *
 * @param Field \Ref{LinBox} field
 */

template <class _Field>
class DenseMatrix : public DenseMatrixBase<typename _Field::Element> 
{
    public:

	typedef typename MatrixCategories::RowMatrixTag<MatrixTraits<DenseMatrix<_Field> > > MatrixCategory;
	typedef _Field Field;
	typedef typename Field::Element   Element;


	DenseMatrix (const Field& F) :  _F(F) , _VD(F) {}

	/** Constructor of a m by n matrix with initial entries which are the 
	 * default constructor value of the field's element type.
	 * @param  F the field of entries; passed so that arithmetic may be done on elements. 
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	DenseMatrix (const Field &F, size_t m, size_t n)
		: DenseMatrixBase<Element> (m, n), _F (F), _VD (F)
	{}

#ifdef __LINBOX_XMLENABLED

	// __LINBOX_XML reader constructor.  Constructs field as well
	DenseMatrix(Reader &R) : DenseMatrixBase<typename Field::Element>(R), _F(R.Down(1)), _VD(_F)
	{ R.Up(1); }


#endif

	/** Constructor of a m by n matrix with entries created by a random iterator.
	 * @param  F the field of entries; passed so that arithmetic may be done on elements. 
	 * @param  m  row dimension
	 * @param  n  column dimension
	 * @para iter, random iterator
	 */
	template<class RandIter>
	DenseMatrix (const Field &F, size_t m, size_t n, const RandIter &iter)
		: DenseMatrixBase<Element> (m, n), _F (F), _VD (F)
	{
		for (typename std::vector<typename Field::Element>::iterator p = _rep.begin (); p != _rep.end (); ++p)
			iter.random (*p);
	}
    
	/** Constructor using a finite vector stream (stream of the rows).
	 * @param  F The field of entries; passed so that arithmetic may be done
	 *           on elements. 
	 * @param  stream A vector stream to use as a source of vectors for this
	 *                matrix
	 */
	template <class StreamVector>
	DenseMatrix (const Field &F, VectorStream<StreamVector> &stream)
		: DenseMatrixBase<Element> (stream.size (), stream.dim ()), _F (F), _VD (F)
	{
		StreamVector tmp;
		typename DenseMatrixBase<Element>::RowIterator p;

		VectorWrapper::ensureDim (tmp, stream.dim ());

		for (p = rowBegin (); p != rowEnd (); ++p) {
			stream >> tmp;
			_VD.copy (*p, tmp);
		}
	}

	/** Constructor from a DenseMatrixBase. Copies all matrix data.
	 * @param F Field over which this matrix' arithmetic will be.
	 * @param M This will contain a complete copy of \Ref{DenseMatrixBase} M.
	 */
	DenseMatrix (const Field &F, DenseMatrixBase<Element> &M)
		: DenseMatrixBase<Element> (M), _F (F), _VD (F)
	{}

	/// Copies {\it all} matrix data.
	DenseMatrix (const DenseMatrix &M)
		: DenseMatrixBase<Element> (M), _F (M._F), _VD (M._F)
	{}

	/** Assignment operator makes a complete copy.
	 */
	DenseMatrix<Field>& operator= (const DenseMatrix<Field>& M) {
		(*this)._rep  = M._rep;
		(*this)._rows = M._rows;
		(*this)._cols = M._cols;
		(*this)._VD   = M._VD;
		const_cast<Field&>((*this)._F)    = M._F;
		return (*this);
	}

	
	/*- Get the number of rows in the matrix
	 * @return Number of rows in matrix
	 */
	size_t rowdim () const
		{ return DenseMatrixBase<Element>::rowdim (); }

	/*- Get the number of columns in the matrix
	 * @return Number of columns in matrix
	 */
	size_t coldim () const
		{ return DenseMatrixBase<Element>::coldim (); }

	/** Retrieve the field over which this matrix is defined
	 * @return Reference to the underlying field
	 */
	const Field &field () const
		{ return _F;}

#ifndef __LINBOX_XMLENABLED

	/*- @name Input and output
	 */

	//@{

	/** Read the matrix from an input stream
	 * @param file Input stream from which to read
	 */
	void read (std::istream &is)
		{ return DenseMatrixBase<Element>::read (is, _F); }
    
	/** Write the matrix to an output stream
	 * @param os Output stream to which to write
	 */
	std::ostream &write (std::ostream &os = std::cout) const
		{ return DenseMatrixBase<Element>::write (os, _F); }
 
	//@}

#else
	ostream &write(ostream &out) const
	{
		Writer W;
		if( toTag(W)) 
			W.write(out);

		return out;
	}

	bool toTag(Writer &W) const
	{
		if(DenseMatrixBase<Element>::toTag(W)) {
			W.insertTagChild();
			_F.toTag(W);
			W.upToParent();
			return true;
		}
		else
			return false;
	}


#endif


	/*- @name Black box interface
	 */

	//@{

	/** Generic matrix-vector apply
	 * y = A * x.
	 * This version of apply allows use of arbitrary input and output vector
	 * types.
	 * @param y Output vector
	 * @param x Input vector
	 * @return Reference to output vector
	 */
	template<class Vect1, class Vect2>
	Vect1 &apply (Vect1 &y, const Vect2 &x) const;

	/** Generic in-place apply
	 * y = A * y.
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

    protected:

	const Field          _F;
	VectorDomain<Field>   _VD;
};

}

#include "dense.inl"

#endif // __BLACKBOX_DENSE_H
