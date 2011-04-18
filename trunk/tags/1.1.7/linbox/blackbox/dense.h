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

#ifndef __LINBOX_blackbox_dense_H
#define __LINBOX_blackbox_dense_H

#include <iostream>
#include <vector>
#include <fstream>

#include "linbox/vector/subiterator.h"
#include "linbox/vector/subvector.h"
#include "linbox/vector/stream.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/matrix/dense.h"
#include <linbox/matrix/matrix-domain.h>
#include <linbox/blackbox/blackbox-interface.h>
#include <linbox/blackbox/factory.h>
#include <linbox/field/hom.h>
#include "linbox/util/matrix-stream.h"

#ifdef __LINBOX_PARALLEL
#include <linbox/blackbox/blackbox_parallel.h>
#endif

namespace LinBox
{

/** \brief Blackbox interface to dense matrix representation. 

 * This is a class of dense matrices
 * templatized by the \link LinBox::FieldArchetype field \endlink 
in which the elements
 * reside. The matrix is stored as a one dimensional STL vector of
 * the elements, in row major order. The interface provides for iteration
 * over rows and over columns.
 *
 * The class conforms to the \link Archetypes archetype \endlink for
\link LinBox::BlackboxArchetype blackboxes \endlink and for the 
\link LinBox::DenseMatrixBase dense matrix container \endlink.
 *
 * Currently, only dense vectors are supported when doing matrix-vector
 * applies.
 *
 * @param _Field a LinBox field class
\ingroup blackbox
 */

template <class _Field>
class DenseMatrix : public BlackboxInterface, public DenseMatrixBase<typename _Field::Element> 
{
    public:

#ifdef __LINBOX_PARALLEL
	BB_list_list sub_list;
#endif

	typedef _Field Field;
	typedef typename Field::Element   Element;
        typedef DenseMatrix<_Field> Self_t;
	typedef typename DenseMatrixBase<typename _Field::Element>::RawIterator RawIterator;
	typedef typename DenseMatrixBase<typename _Field::Element>::ConstRawIterator ConstRawIterator;
	typedef typename DenseMatrixBase<typename _Field::Element>::RawIndexedIterator RawIndexedIterator;
	typedef typename DenseMatrixBase<typename _Field::Element>::ConstRawIndexedIterator ConstRawIndexedIterator;

	DenseMatrix (const Field& F) :  _F(F) , _MD(F), _AT (*this) {}

	/** Constructor of a m by n matrix with initial entries which are the 
	 * default constructor value of the field's element type.
	 * @param  F the field of entries; passed so that arithmetic may be done on elements. 
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	DenseMatrix (const Field &F, size_t m, size_t n)
		: DenseMatrixBase<Element> (m, n), _F (F), _MD (F), _AT (*this)
	{}

	/** Constructor of a m by n matrix with entries created by a random iterator.
	 * @param  F the field of entries; passed so that arithmetic may be done on elements. 
	 * @param  m  row dimension
	 * @param  n  column dimension
	 * @para iter, random iterator
	 */
	template<class RandIter>
	DenseMatrix (const Field &F, size_t m, size_t n, const RandIter &iter)
		: DenseMatrixBase<Element> (m, n), _F (F), _MD (F), _AT (*this)
	{
		for (typename std::vector<typename Field::Element>::iterator p = DenseMatrixBase<Element>::_rep.begin (); p != DenseMatrixBase<Element>::_rep.end (); ++p)
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
		: DenseMatrixBase<Element> (stream.size (), stream.dim ()), _F (F), _MD (F), _AT (*this)
	{
		StreamVector tmp;
		typename DenseMatrixBase<Element>::RowIterator p;

		VectorWrapper::ensureDim (tmp, stream.dim ());

		VectorDomain<Field> _VD(F);

		for (p = DenseMatrixBase<Element>::rowBegin (); p != DenseMatrixBase<Element>::rowEnd (); ++p) {
			stream >> tmp;
			_VD.copy (*p, tmp);
		}
	}

	/** Constructor from a MatrixStream
	 * @param ms A matrix stream properly initialized
	 */
	DenseMatrix( MatrixStream<Field>& ms )
		:DenseMatrixBase<Element>(ms), _F( ms.getField() ), _MD( ms.getField() ), _AT(*this)
	{ }

	/** Constructor from a DenseMatrixBase. Copies all matrix data.
	 * @param F Field over which this matrix' arithmetic will be.
	 * @param M This will contain a complete copy of \ref{DenseMatrixBase} M.
	 */
	DenseMatrix (const Field &F, DenseMatrixBase<Element> &M)
		: DenseMatrixBase<Element> (M), _F (F), _MD (F), _AT (*this)
	{}

	/// Copies {\it all} matrix data.
	DenseMatrix (const DenseMatrix &M)
		: DenseMatrixBase<Element> (M), _F (M._F), _MD (M._F), _AT (*this)
	{}

    
	/** Assignment operator makes a complete copy.
	 */
	DenseMatrix<Field>& operator= (const DenseMatrix<Field>& M) {
		(*this)._rep  = M._rep;
		(*this)._rows = M._rows;
		(*this)._cols = M._cols;
		(*this)._MD   = const_cast<MatrixDomain<Field>&>(M._MD);
		(*this)._F    = const_cast<Field&>(M._F);
		return (*this);
	}


	template<typename _Tp1>
	struct rebind
	{ 
		typedef DenseMatrix<_Tp1> other; 
		
		void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
// 			Ap = new other(F, A.rowdim(), A.coldim());
			typename Self_t::ConstRawIterator A_p;
			typename other::RawIterator Ap_p;
			Hom<Field, _Tp1> hom(A. field(), F);
			for (A_p = A. rawBegin(), Ap_p = Ap.rawBegin();
			     A_p != A. rawEnd(); ++ A_p, ++ Ap_p) 
				hom.image (*Ap_p, *A_p);
		}
	};

    	template<typename _Tp1>
	DenseMatrix (const DenseMatrix<_Tp1> &M, const Field& F)
		: DenseMatrixBase<Element> (M.rowdim(),M.coldim()), _F (F), _MD (F), _AT (*this) {
            typename DenseMatrix<_Tp1>::template rebind<Field>()(*this, M, F);
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

	Field &field () 
		{ return _F;}

	/*- @name Input and output
	 */

	//@{

	/** Read the matrix from an input stream
	 * @param file Input stream from which to read
	 */
	std::istream& read (std::istream &is)
	{ return DenseMatrixBase<Element>::read (is, _F); }
    
	/** Write the matrix to an output stream
	 * @param os Output stream to which to write
	 */
	std::ostream &write (std::ostream &os = std::cout) const
		{ return DenseMatrixBase<Element>::write (os, _F); }
 
	//@}


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
  
    
	// destructor
	~DenseMatrix ( ) {

#ifdef __LINBOX_PARALLEL

                BB_list_list::iterator p;

                BB_list::iterator e_p;

                for (p = sub_list. begin(); p != sub_list. end(); ++ p)

                        for (e_p = p -> second. begin();
                             e_p != p -> second. end(); ++ e_p) {

                                Thread::terminate_thread (*e_p);

                                delete (*e_p);
                          }

#endif
	}

	//@}

    protected:

	//const Field          _F;
	Field          _F;
	MatrixDomain<Field>   _MD;
	TransposeMatrix<DenseMatrix<Field> > _AT;
};

template <class Field>
struct MatrixTraits< DenseMatrix<Field> >
{
	typedef DenseMatrix<Field> MatrixType;
	typedef typename MatrixCategories::RowMatrixTag MatrixCategory;
};

template <class Field>
struct MatrixTraits< const DenseMatrix<Field> >
{
	typedef const DenseMatrix<Field> MatrixType;
	typedef typename MatrixCategories::RowMatrixTag MatrixCategory;
};

/** Dense matrix factory
  * This class ingerits \ref{BlackboxFactory} and provides a method for using a
  * \ref{DenseMatrixBase} object with integer or rational data type as input to
  * the high-level intger and rational solutions functions.
  */

template< class Field, class BElement >
class DenseMatrixFactory : public BlackboxFactory<Field,DenseMatrix<Field> >
{
    private:
    	const DenseMatrixBase<BElement>& _A;

    public:
    	DenseMatrixFactory( const DenseMatrixBase<BElement> &A ) :_A(A) {}

	DenseMatrix<Field>* makeBlackbox( const Field& F );

	integer& maxNorm( integer& res );

	integer& hadamardBound( integer& res ) const;

	size_t rowdim() { return _A.rowdim(); }
	size_t coldim() { return _A.coldim(); }
};

}

#include "dense.inl"

#endif // __LINBOX_blackbox_dense_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
