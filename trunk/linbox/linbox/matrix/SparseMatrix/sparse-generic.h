/* linbox/matrix/sparse.h
 * Copyright (C) 2001-2002 Bradford Hovinen
 *               1999-2001 William J Turner,
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 * 2003-01-11  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Move from blackbox/sparse-base.h to matrix/sparse.h
 * ------------------------------------
 * 2002-11-28  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 *   - Renamed ColOfRowsIterator to RowIterator
 *   - Named template argument _Row rather than Row; add a typedef to Row
 *   - Named template argument _Element rather than Row; add a typedef to Element
 *   - Renamed IndexIterator as IndexedIterator, and adjusted to match
 *     interface in DenseMatrixBase
 * ------------------------------------
 * 2002-08-06  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Renamed to sparse-base.h from sparse0-base.h
 * ------------------------------------
 * Modified by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Refactoring:
 *   - Eliminated SparseMatrixAux and moved that functionality into Sparse0
 *   - Made SparseMatrix parameterized only on the element type
 *   - New read/write implementations for SparseMatrix, supporting multiple
 *     formats
 *   - Eliminated Gaussian elimination code
 *   - Added iterators, including ColOfRowsIterator, Iterator, and
 *     IndexIterator
 *   - Eliminated operator []; added getEntry; changed put_value to setEntry
 * ------------------------------------
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_matrix_sparse_H
#define __LINBOX_matrix_sparse_H

#ifndef _SP_BB_VECTOR_
#include <vector>
#define _SP_BB_VECTOR_ std::vector
#endif

#include <utility>
#include <iostream>
#include <algorithm>

#include "linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/linbox-tags.h"
#include "linbox/matrix/sparse-formats.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/transpose-matrix.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/solutions/solution-tags.h"
#include "linbox/matrix/matrix-traits.h"
#include "linbox/field/hom.h"

namespace LinBox
{
	template<class Field>
	class MatrixDomain ;

	// Forward declaration
	template <class _Field,
		 class _Row   = typename RawVector<typename _Field::Element>::Sparse,
		 class Trait  = typename VectorTraits<_Row>::VectorCategory>
		 class SparseMatrix;


	// Small helper classes to make read and write easier
	template <class _Field, class Row,
		 class Trait = typename VectorTraits<Row>::VectorCategory>
		 class SparseMatrixWriteHelper {
		 public:
			 typedef         _Field          Field;
			 typedef typename Field::Element Element;

			 static std::ostream &write (const SparseMatrix<Field, Row> &A
						     , std::ostream &os
						     , LINBOX_enum(Tag::FileFormat) format);
		 };

	template <class Field, class Row, class Trait = typename VectorTraits<Row>::VectorCategory>
	class SparseMatrixReadWriteHelper : public SparseMatrixWriteHelper<Field, Row, Trait> {
		static std::istream &readTurner    (SparseMatrix<Field, Row> &A
						    , std::istream &is
						    // , const Field &F
						    , char *buf);
		static std::istream &readGuillaume (SparseMatrix<Field, Row> &A
						    , std::istream &is
						    // , const Field &F
						    , char *buf);
		static std::istream &readMatlab    (SparseMatrix<Field, Row> &A
						    , std::istream &is
						    // , const Field &F
						    , char *buf);
		static std::istream &readPretty    (SparseMatrix<Field, Row> &A
						    , std::istream &is
						    // , const Field &F
						    , char *buf);
		static std::istream &readMagmaCpt  (SparseMatrix<Field, Row> &A
						    , std::istream &is
						    // , const Field &F
						    , char *buf);

	public:
		typedef typename Field::Element Element;

		static std::istream &read (SparseMatrix<Field, Row> &A
					   , std::istream &is
					   , LINBOX_enum(Tag::FileFormat) format);
	};

	/** Sparse matrix container
	 * This class acts as a generic row-wise container for sparse
	 * matrices. It is designed to provide various methods to access the
	 * entries of the matrix.
	 * It meets the blackbox archetype
	 * @tparam _Field Field type
	 * @tparam _Row   LinBox sparse vector type to use for rows of matrix
	 \ingroup matrix
	 */
	template <class _Field, class _Row, class Trait>
	class SparseMatrix {
	public:

		typedef         _Field          Field;
		typedef typename Field::Element Element;
		typedef _Row Row;
		typedef const Row ConstRow;
		typedef typename _SP_BB_VECTOR_<Row> Rep;
		typedef SparseMatrix<_Field, _Row, Trait> Self_t;

#ifdef __LINBOX_PARALLEL
		BB_list_list sub_list;
#endif



		/** Constructor.
		 * Note: the copy constructor and operator= will work as intended
		 *       because of STL's container design
		 * @param  m  row dimension
		 * @param  n  column dimension
		 */
		SparseMatrix (const Field & F,size_t m, size_t n) :
			_field (F),
			_MD(F),_AT(*this),
			_matA(m), _m(m), _n(n)
		{};

		SparseMatrix (const Field & F) :
			_field (F),
			_MD(F),_AT(*this),
			_matA(0), _m(0), _n(0)
		{};



		/** Constructor from a MatrixStream
		*/
		SparseMatrix ( MatrixStream<Field>& ms );


		/** Copy constructor.
		*/
		SparseMatrix (const SparseMatrix<Field, Row, Trait> &A);

		/** Convert constructor.
		*/
		template<class VectorType>
		SparseMatrix (const SparseMatrix<Field, VectorType, Trait> &A);

		template<typename _Tp1, typename _R1 = typename Rebind<_Row,_Tp1>::other >
		struct rebind {
			typedef SparseMatrix<_Tp1, _R1> other;

			void operator() (other & Ap, const Self_t& A) {

				typename _Tp1::Element e;

				Hom<typename Self_t::Field, _Tp1> hom(A.field(), Ap.field());
				for( typename Self_t::ConstIndexedIterator
				     indices = A.IndexedBegin();
				     (indices != A.IndexedEnd()) ;
				     ++indices ) {
					// hom. image (e, A.getEntry(indices.rowIndex(),indices.colIndex()) );
					hom. image (e, indices.value() );
					if (!Ap.field().isZero(e))
						Ap.setEntry (indices.rowIndex(),
							     indices.colIndex(), e);
				}
			}

		};

		template<typename _Tp1, typename _Rw1, typename _myT>
		SparseMatrix (const SparseMatrix<_Tp1, _Rw1, _myT> &Mat, const Field& F) :
			_field (&F),
			_MD (F), _AT (*this)
			, _matA (Mat.rowdim()), _m (Mat.rowdim()), _n (Mat.coldim())
		{
				typename SparseMatrix<_Tp1,_Rw1, _myT>::template rebind<Field,_Row>()(*this, Mat);
		}

		template<class VectStream>
		SparseMatrix (const Field &F, VectStream &stream) :
			_field (F), _MD (F), _AT (*this)
			, _matA (stream.size()), _m (stream.size()), _n (stream.dim())
		{
			typename Self_t::RowIterator i;

			for (i = Self_t::rowBegin (); i != Self_t::rowEnd (); ++i)
				stream >> *i;
		}



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

		/** Retreive row dimension of the matrix.
		 * @return integer number of rows of SparseMatrix matrix.
		 */
		size_t rowdim () const
		{
			return _m;
		}

		/** Retreive column dimension of matrix.
		 * @return integer number of columns of SparseMatrix matrix.
		 */
		size_t coldim () const
		{
			return _n;
		}

		/** Retreive number of elements in the matrix.
		 * @return integer number of elements of SparseMatrix matrix.
		 * @bug should it be elements or non zero elements ? @see ELL
		 */
		size_t size () const
		{
			size_t s(0);
			for(typename Rep::const_iterator it = _matA.begin(); it != _matA.end(); ++it)
				s+= LinBox::RawVector<Element>::size(*it);
			return s;
		}

		/** Read a matrix from the given input stream using field read/write
		 * @param is Input stream from which to read the matrix
		 * @param format Format of input matrix
		 */
		std::istream &read (std::istream &is,   LINBOX_enum(Tag::FileFormat) format /*= Tag::FileFormat::Detect*/);

		// Read from matrix market format
		std::istream &read (std::istream &is)
		{
			MatrixStream<Field> ms(field(), is);
			if( !ms.getDimensions( this->_m, this->_n ) )
				throw ms.reportError(__func__,__LINE__);
			this->_matA.resize( this->_m );
			Element val;
			size_t i, j;
			while( ms.nextTriple(i,j,val) ) {
				setEntry(i,j,val);
			}
			if( ms.getError() > END_OF_MATRIX )
				throw ms.reportError(__func__,__LINE__);
			return is;
		}


		/// Write in matrix market format
		std::ostream &write (std::ostream &os) const
		{
			// typedef SparseMatrixBase<Element, _Row> SMB;
			writeMMCoordHeader(os, *this, this->size(), "SparseMatrix");
			return this->write(os, Tag::FileFormat::OneBased);
		}

		/** Write a matrix to the given output stream using field read/write
		 * @param os Output stream to which to write the matrix
		 * @param F Field with which to write
		 * @param format Format with which to write
		 */
		std::ostream &write (std::ostream &os,  LINBOX_enum(Tag::FileFormat) format /*  = Tag::FileFormat::Pretty*/) const;

		/** Set an individual entry
		 * Setting the entry to 0 will remove it from the matrix
		 * @param i Row index of entry
		 * @param j Column index of entry
		 * @param value Value of the new entry
		 */
		void setEntry (size_t i, size_t j, const Element &value);
		void finalize(){}

		/** Get a writeable reference to an entry in the matrix
		 * If there is no entry at the position (i, j), then a new entry
		 * with a value of zero is inserted and a reference  to it is
		 * returned.
		 * @param i Row index of entry
		 * @param j Column index of entry
		 * @return Reference to matrix entry
		 */
		Element &refEntry (size_t i, size_t j);

		/** Get a read-only individual entry from the matrix
		 * @param i Row index
		 * @param j Column index
		 * @return Const reference to matrix entry
		 */
		const Element &getEntry (size_t i, size_t j) const;

		/** Get an entry and store it in the given value
		 * This form is more in the Linbox style and is provided for interface
		 * compatibility with other parts of the library
		 * @param x Element in which to store result
		 * @param i Row index
		 * @param j Column index
		 * @return Reference to x
		 */
		Element &getEntry (Element &x, size_t i, size_t j) const;

		/** @name Columns of rows iterator
		 * The columns of row iterator gives each of the rows of the
		 * matrix in ascending order. Dereferencing the iterator yields
		 * a row vector in sparse sequence format
		 */

		typedef typename Rep::iterator RowIterator;
		typedef typename Rep::const_iterator ConstRowIterator;

		RowIterator rowBegin ();
		RowIterator rowEnd ();
		ConstRowIterator rowBegin () const;
		ConstRowIterator rowEnd () const;

		/** @name Raw iterator
		 * The raw iterator is a method for accessing all nonzero
		 * entries in the matrix in some unspecified order. This can be
		 * used, e.g. to reduce all matrix entries modulo a prime before
		 * passing the matrix into an algorithm.
		 */

		class Iterator;
		class ConstIterator;

		/// Begin.
		Iterator Begin ();
		/// End.
		Iterator End ();
		/// const Begin.
		ConstIterator Begin () const;
		/// const End
		ConstIterator End () const;

		/** @name Index iterator
		 * The index iterator gives the row, column indices of all matrix
		 * elements in the same order as the raw iterator above. Its value type
		 * is an STL pair with the row and column indices, starting at 0, in the
		 * first and second positions, respectively.
		 */

		class IndexedIterator;
		class ConstIndexedIterator;

		/// IndexedBegin
		IndexedIterator IndexedBegin ();
		/// IndexedEnd
		IndexedIterator IndexedEnd ();
		/// const IndexedBegin
		ConstIndexedIterator IndexedBegin () const;
		/// const IndexedEnd
		ConstIndexedIterator IndexedEnd () const;

		/** Retrieve a row as a writeable reference
		 * @param i Row index
		 */
		Row &getRow (size_t i);

		/** Retrieve a row as a writeable reference.
		 * @param i Row index
		 */
		Row &operator [] (size_t i);

		/** Retrieve a row as a read-only reference.
		 * @param i Row index
		 */
		ConstRow &operator [] (size_t i) const;

		/** Compute the column density, i.e the number of entries per column.
		 * @param v Vector in which to store column density
		 */
		template <class Vector>
		Vector &columnDensity (Vector &v) const;

		/** Construct the transpose of this matrix and place it in the
		 * matrix given.
		 * @param AT
		 */
		SparseMatrix &transpose (SparseMatrix &AT) const;

		const Field & field() const
		{
			return _field;
		}

		/** Matrix-vector product
		 * \f$y = A x\f$.
		 * @return reference to output vector y
		 * @param  x input vector
		 * @param y
		 */
		template <class OutVector, class InVector>
		OutVector &apply (OutVector &y, const InVector &x) const
		{
#ifdef __LINBOX_PARALLEL
			return BlackboxParallel (y, *this, x, BBBase::Apply);
#else
			return _MD.vectorMul (y, *this, x);
#endif
		}

		/** Transpose matrix-vector product
		 * \f$ y = A^T x\f$.
		 * @return reference to output vector y
		 * @param  x input vector
		 * @param y
		 */
		template <class OutVector, class InVector>
		OutVector &applyTranspose (OutVector& y, const InVector &x) const
		{
#ifdef __LINBOX_PARALLEL
			return BlackboxParallel (y, *this, x, BBBase::ApplyTranspose);
#else
			return _MD.vectorMul (y, _AT, x);
#endif
		}

		const Rep & getRep() const
		{
			return _matA;
		}

	protected:

		friend class SparseMatrixWriteHelper<Field, Row>;
		friend class SparseMatrixReadWriteHelper<Field, Row>;

		Rep               _matA;
		size_t            _m;
		size_t            _n;
		const Field & _field;

		MatrixDomain<Field>       _MD; // Matrix domain for matrix operations
		TransposeMatrix<SparseMatrix<_Field, _Row> > _AT;


		// template<class F, class R, class T> friend class SparseMatrix;
	};

} // LinBox

#include "linbox/matrix/SparseMatrix/sparse-sequence-vector.h"
#include "linbox/matrix/SparseMatrix/sparse-parallel-vector.h"
#include "linbox/matrix/SparseMatrix/sparse-associative-vector.h"

#include "linbox/matrix/SparseMatrix/sparse-generic.inl"

namespace LinBox {

	template <class Field, class Row>
	std::ostream &operator << (std::ostream &os, const SparseMatrix<Field, Row> &A)
	{
		return A.write (os);
	}

	template <class Field, class Row>
	std::istream &operator >> (std::istream &is, SparseMatrix<Field, Row> &A)
	{
		return A.read (is);
	}

	template <class Field, class Row, class Trait>
	struct MatrixTraits< SparseMatrix<Field, Row, Trait> >
	{
		typedef SparseMatrix<Field, Row, Trait> MatrixType;
		typedef typename MatrixCategories::RowMatrixTag MatrixCategory;
	};

	template <class Field, class Row, class Trait>
	struct MatrixTraits< const SparseMatrix<Field, Row, Trait> >
	{
		typedef const SparseMatrix<Field, Row, Trait> MatrixType;
		typedef typename MatrixCategories::RowMatrixTag MatrixCategory;
	};

	template<class A, class B, class C>
	struct GetEntryCategory<SparseMatrix<A,B,C> > {
		typedef SolutionTags::Local Tag;
	} ;

	template <class Field, class _Row>
	struct MatrixTraits< SparseMatrix<Field, _Row> >
	{
		typedef SparseMatrix<Field, _Row> MatrixType;
		typedef MatrixCategories::RowMatrixTag MatrixCategory;
	};

	template<class A, class B>
	struct GetEntryCategory<SparseMatrix<A,B> >
	{
		typedef SolutionTags::Local Tag;
	};

} // namespace LinBox


#endif // __LINBOX_matrix_sparse_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
