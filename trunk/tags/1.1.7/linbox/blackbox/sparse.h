/* linbox/blackbox/sparse.h
 * Copyright (C) 1999-2005 William J Turner,
 *               2001-2002 Bradford Hovinen
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-08-06  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Renamed to sparse.h from sparse0.h
 * ------------------------------------
 * Modified by Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 * 	       28.08.2002 : added back : field()
 *
 * Refactoring:
 *   - Eliminated SparseMatrixAux and moved that functionality into Sparse0
 *   - Made SparseMatrixBase parameterized only on the element type
 *   - New read/write implementations for SparseMatrixBase, supporting multiple
 *     formats
 *   - Eliminated Gaussian elimination code
 *   - Added iterators, including ColOfRowsIterator, RawIterator, and
 *     RawIndexIterator
 *   - Eliminated operator []; added getEntry; changed put_value to setEntry
 * ------------------------------------
 * Modified by W. J. Turner <wjturner@acm.org>
 *	24.06.2005 : Removed using declarations
 * ------------------------------------
 * 
 * See COPYING for license information.
 */

#ifndef __LINBOX_blackbox_sparse_H
#define __LINBOX_blackbox_sparse_H

#include "linbox/linbox-config.h"
#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/matrix/sparse.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/vector/stream.h"
#include "linbox/util/field-axpy.h"
#include <linbox/field/hom.h>
#include <linbox/field/rebind.h>

namespace LinBox
{

/** \brief vector of sparse rows.
    
 * This is a generic black box for a sparse matrix. It inherits
 * LinBox::SparseMatrixBase, which implements all of the underlying
 * accessors and iterators.
 * \ingroup blackbox
 */
template <class _Field,
	  class _Row    = typename LinBox::Vector<_Field>::Sparse>
class SparseMatrix : public BlackboxInterface, public SparseMatrixBase<typename _Field::Element, _Row> 
{
    public:

	typedef _Field Field;
	typedef typename Field::Element Element;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::Row Row;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::Rep Rep;
        typedef SparseMatrix<_Field, _Row> Self_t;
    

#ifdef __LINBOX_PARALLEL
	BB_list_list sub_list;
#endif

	FileFormatTag Format;

	typedef typename SparseMatrixBase<typename Field::Element, _Row>::RawIterator RawIterator;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::RawIndexedIterator RawIndexedIterator;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::ConstRawIterator ConstRawIterator;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::ConstRawIndexedIterator ConstRawIndexedIterator;

	/** Constructor.
	 * Builds a zero m x n matrix
	 * Note: the copy constructor and operator= will work as intended
	 *       because of STL's container design
	 * @param  F  Field over which entries exist
	 * @param  m  Row dimension
	 * @param  n  Column dimension
	 */
// 	SparseMatrix (const Field &F)
// 		: SparseMatrixBase<Element, _Row> (0,0), _F (F), _VD (F), _MD (F), _AT (*this)
// 	{            std::cerr << "default cstor" << std::endl;
// }

// 	SparseMatrix (const Field &F, size_t m, size_t n)
// 		: SparseMatrixBase<Element, _Row> (m, n), _F (F), _VD (F), _MD (F), _AT (*this)
// 	{            std::cerr << "default cstor : " <<  m << "x" << n << std::endl;
// }
	SparseMatrix (const Field &F, size_t m=0, size_t n=0)
		: SparseMatrixBase<Element, _Row> (m, n), _F (F), _VD (F), _MD (F), _AT (*this)
	{ }

	/** Constructor from a vector stream
	 * @param  F  Field over which entries exist
	 * @param  stream  Stream with which to generate row vectors
	 */
        template<class VectStream>
	SparseMatrix (const Field &F, VectStream &stream)
		: SparseMatrixBase<Element, _Row> (stream.size (), stream.dim ()),
		  _F (F), _VD (F), _MD (F), _AT (*this)
	{
		typename SparseMatrixBase<Element, _Row>::RowIterator i;

		for (i = SparseMatrixBase<Element, _Row>::rowBegin (); i != SparseMatrixBase<Element, _Row>::rowEnd (); ++i)
			stream >> *i;
	}


	/** Constructor from a MatrixStream
	 * @param ms A matrix stream properly initialized
	 */
	SparseMatrix( MatrixStream<Field>& ms )
		:SparseMatrixBase<Element,_Row>(ms), _F(ms.getField()), _VD(ms.getField()), _MD(ms.getField()), _AT(*this)
	{ }

	/** Copy constructor
	 */
	SparseMatrix (const SparseMatrix<Field, Row> &B)
		: SparseMatrixBase<Element, _Row> (B), _F (B._F), _VD (B._F), _MD (B._F), _AT (*this)
	{ }

	/** Row type Converter constructor
	 */
        template<class VectorType>
	SparseMatrix (const SparseMatrix<Field, VectorType> &B)
		: SparseMatrixBase<Element, _Row> (B), _F (B._F), _VD (B._F), _MD (B._F), _AT (*this)
	{ }

	/** Destructor. */
	~SparseMatrix () {
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

	/** Matrix-vector product
	 * y = A x.
	 * @return reference to output vector y
	 * @param  x input vector
	 */
	template <class OutVector, class InVector>
	OutVector &apply (OutVector &y, const InVector &x) const {
#ifdef __LINBOX_PARALLEL
		return BlackboxParallel (y, *this, x, BBBase::Apply);
#else
		return _MD.vectorMul (y, *this, x);
#endif
	}

	/** Transpose matrix-vector product
	 * y = A^T x.
	 * @return reference to output vector y
	 * @param  x input vector
	 */
	template <class OutVector, class InVector>
	OutVector &applyTranspose (OutVector& y, const InVector &x) const { 
#ifdef __LINBOX_PARALLEL
		return BlackboxParallel (y, *this, x, BBBase::ApplyTranspose);
#else
		return _MD.vectorMul (y, _AT, x);
#endif
	}


	template<typename _Tp1, typename _Rw1 = typename Rebind<_Row, _Tp1>::other> 
	struct rebind { 
		typedef SparseMatrix<_Tp1, _Rw1> other;	
		
		void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
// 			Ap = new other(F, A.rowdim(), A.coldim());
	
			typename _Tp1::Element e;
			Hom<typename Self_t::Field, _Tp1> hom(A.field(), F);
			for( typename Self_t::ConstRawIndexedIterator
                                 indices = A.rawIndexedBegin();
                             	 (indices != A.rawIndexedEnd()) ; 
                                 ++indices ) {
//                             hom. image (e, A.getEntry(indices.rowIndex(),indices.colIndex()) );
                            hom. image (e, indices.value() );
                            if (!F.isZero(e)) 
                                Ap.setEntry (indices.rowIndex(), 
                                                indices.colIndex(), e);
			}
		}
	};

    	template<typename _Tp1, typename _Rw1>
	SparseMatrix (const SparseMatrix<_Tp1, _Rw1> &M, const Field& F)
		: SparseMatrixBase<Element, _Row> (M.rowdim(),M.coldim()), _F (F), _VD (F), _MD (F), _AT (*this) {
            typename SparseMatrix<_Tp1,_Rw1>::template rebind<Field,_Row>()(*this, M, F);
        }




	/** Retreive row dimensions of Sparsemat matrix.
	 * @return integer number of rows of SparseMatrix0Base matrix.
	 */
	size_t rowdim () const { return SparseMatrixBase<Element, _Row>::_m; }

	/** Retreive column dimensions of Sparsemat matrix.
	 * @return integer number of columns of SparseMatrix0Base matrix.
	 */
	size_t coldim () const { return SparseMatrixBase<Element, _Row>::_n; }

	/** Read the matrix from a stream in the given format
	 * @param is Input stream from which to read the matrix
	 * @param format Format of input matrix
	 * @return Reference to input stream
	 */
	std::istream &read (std::istream &is, FileFormatTag format = FORMAT_DETECT)
	{ return SparseMatrixBase<Element, _Row>::read (is, _F, format); }

	/** Write the matrix to a stream in the given format
	 * @param os Output stream to which to write the matrix
	 * @param format Format of output
	 * @return Reference to output stream
	 */
	std::ostream &write (std::ostream &os, FileFormatTag format = FORMAT_PRETTY) const
		{ return SparseMatrixBase<Element, _Row>::write (os, _F, format); }

	// JGD 28.08.2002
	/** Access to the base field
	 */
	const Field& field () const { return _F;}

    protected:

	const Field                             _F;      // Field used for all arithmetic
	VectorDomain<Field>                     _VD;     // Vector domain for matrix operations
	MatrixDomain<Field>                     _MD;     // Matrix domain for matrix operations

	TransposeMatrix<SparseMatrix<_Field, _Row> > _AT;

    	template<class F, class R> friend class SparseMatrix;
};

/** Sparse matrix factory
 * This class inherits \ref{BlackboxFactory} and provides a method for using a
 * \ref{SparseMatrixBase} object with integer or rational data type as input to
 * the high-level integer and rational solutions functions.
 */

template <class Field,
	  class BElement = typename Field::Element,
	  class Row      = typename LinBox::Vector<Field>::Sparse,
	  class BRow     = typename LinBox::RawVector<BElement>::Sparse>
class SparseMatrixFactory : public BlackboxFactory<Field,SparseMatrix<Field,Row> > 
{//otot
	const SparseMatrixBase<BElement, BRow> &_A;

    public:

	SparseMatrixFactory (const SparseMatrixBase<BElement, BRow> &A)
		: _A (A) 
	{}

	// FIXME: This function assumes basically that the matrix is over the integers

	SparseMatrix<Field,Row> *makeBlackbox (const Field &F);

	integer &maxNorm (integer &res)
	{
		typename SparseMatrixBase<BElement, BRow>::ConstRawIterator i;

		res = 0L;

		integer tmp;

		for (i = _A.rawBegin (); i != _A.rawEnd (); ++i) {
			tmp = abs (*i);

			if (res < tmp)
				res = tmp;
		}

		return res;
	}

	size_t rowdim ()
        	{ return _A.rowdim (); }
	size_t coldim ()
		{ return _A.coldim (); }

	// A better bound for determinant of an integer sparse matrix, ZW
	integer &hadamardBound (integer& res) const {
		return hadamardBound (res, typename VectorTraits<Row>::VectorCategory());
	}

	integer &hadamardBound (integer& res,  VectorCategories::SparseParallelVectorTag) const {
		typedef typename SparseMatrixBase<BElement, BRow>::ConstRowIterator RowIterator;
		typedef typename SparseMatrixBase<BElement, BRow>::ConstRow::second_type::const_iterator EltIterator;

		res = 1L;

		integer tmp;
		RowIterator row_p;
		EltIterator elt_p;

		for (row_p = _A. rowBegin(); row_p != _A. rowEnd(); ++ row_p) {
			tmp = 0;

			for (elt_p = row_p -> second. begin(); elt_p != row_p -> second. end(); ++ elt_p)
				tmp += static_cast<integer>(*elt_p) * (*elt_p);

			res *=tmp;
		}

		res = sqrt (res);
		return res;
	}

	integer &hadamardBound (integer& res,  VectorCategories::SparseSequenceVectorTag) const{
		typedef typename SparseMatrixBase<BElement, BRow>::ConstRowIterator RowIterator;
		typedef typename SparseMatrixBase<BElement, BRow>::ConstRow::const_iterator EltIterator;

		res = 1L;
		 
		integer tmp;
		RowIterator row_p;
		EltIterator elt_p;
		 
		for ( row_p = _A. rowBegin(); row_p != _A. rowEnd(); ++ row_p) {
			tmp = 0;
		      
			for (EltIterator elt_p = row_p -> begin(); elt_p != row_p -> end(); ++ elt_p)
				tmp += static_cast<integer>(elt_p -> second) * (elt_p -> second);
		
			res *=tmp; 
		} 
		 
		res = sqrt (res); 
		return res; 
	} 

	integer &hadamardBound (integer& res,  VectorCategories::SparseAssociativeVectorTag) const {
		typedef typename SparseMatrixBase<BElement, BRow>::ConstRowIterator RowIterator;
		typedef typename SparseMatrixBase<BElement, BRow>::ConstRow::const_iterator EltIterator;

		res = 1L;
		 
		integer tmp;
		RowIterator row_p;
		EltIterator elt_p;
		 
		for ( row_p = _A. rowBegin(); row_p != _A. rowEnd(); ++ row_p) {
			tmp = 0;

			for (elt_p = row_p -> begin(); elt_p != row_p -> end(); ++ elt_p)
				tmp += static_cast<integer>(elt_p -> second) * (elt_p -> second);
			 
			res *=tmp; 
		} 
		 
		res = sqrt (res); 
		return res; 
	} 
};

#if !defined(__INTEL_COMPILER) && !defined(__CUDACC__)
template <>
#endif
template <class Field, class _Row>
struct MatrixTraits< SparseMatrix<Field, _Row> >
{ 
	typedef SparseMatrix<Field, _Row> MatrixType;
	typedef MatrixCategories::RowMatrixTag MatrixCategory;
};

template <class Field, class _Row>
struct MatrixTraits< const SparseMatrix<Field, _Row> >
{ 
	typedef const SparseMatrix<Field, _Row> MatrixType;
	typedef MatrixCategories::RowMatrixTag MatrixCategory;
};

} // namespace LinBox

#include "linbox/blackbox/sparse.inl"

#endif // __LINBOX_blackbox_sparse_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
