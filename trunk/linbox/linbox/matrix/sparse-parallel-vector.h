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

#ifndef __LINBOX_matrix_sparse_parallel_H
#define __LINBOX_matrix_sparse_parallel_H

namespace LinBox
{


	// Specialization of the above for sparse parallel vectors
	template <class _Field, class Row>
	class SparseMatrixWriteHelper<_Field, Row, VectorCategories::SparseParallelVectorTag > {
	public:
		typedef         _Field          Field;
		typedef typename Field::Element Element;

		static std::ostream &write (const SparseMatrix<Field, Row> &A, std::ostream &os
					    // , const Field &F
					    , LINBOX_enum(Tag::FileFormat) format);
	};


	/* Specialization for sparse parallel vectors */

	template <class _Field, class _Row>
	class SparseMatrix<_Field, _Row, VectorCategories::SparseParallelVectorTag > {
	public:

		typedef         _Field          Field;
		typedef typename Field::Element Element;

		typedef _Row                Row;
		typedef const Row      ConstRow;
		typedef _SP_BB_VECTOR_<Row> Rep;
		typedef VectorCategories::SparseParallelVectorTag myTrait;
		typedef SparseMatrix<_Field, _Row, myTrait> Self_t;

#ifdef __LINBOX_PARALLEL
		BB_list_list sub_list;
#endif


	template<typename _Tp1, typename _R1 = typename Rebind<_Row,_Tp1>::other>
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
	template<typename _Tp1, typename _Rw1>
	SparseMatrix (const SparseMatrix<_Tp1, _Rw1, myTrait> &Mat, const Field& F) :
		_field (F),
		_MD (F), _AT (*this)
			, _matA (Mat.rowdim()), _m (Mat.rowdim()), _n (Mat.coldim())
	{
		typename SparseMatrix<_Tp1,_Rw1, myTrait>::template rebind<Field,_Row>()(*this, Mat);
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

		SparseMatrix (const Field & F, size_t m, size_t n) :
			_field(F),
			_MD(F),_AT(*this),
			_matA (m), _m (m), _n (n)
		{}

	SparseMatrix (const Field & F) :
			_field (F),
			_MD(F),_AT(*this),
			_matA(0), _m(0), _n(0)
		{};

		SparseMatrix (const SparseMatrix<Field, Row> &A) :
			_field(A.field()),
			_MD(A.field()),_AT(*this),
			_matA (A._matA), _m (A._m), _n (A._n)
		{}

		template<class VectorType>
		SparseMatrix (const SparseMatrix<Field, VectorType> &A) :
			_field(A.field()),
			_MD(A.field()),_AT(*this),
			_matA(A._m), _m (A._m), _n (A._n)
		{
			typename Rep::iterator meit = this->_matA.begin();
			typename SparseMatrix<Field, VectorType>::Rep::const_iterator copit = A._matA.begin();
			for( ; meit != this->_matA.end(); ++meit, ++copit)
				LinBox::RawVector<Element>::convert(*meit, *copit);
		}

		/** Constructor from a MatrixStream
		*/
		SparseMatrix ( MatrixStream<Field>& ms );

		~SparseMatrix () {}

		size_t rowdim () const
		{
			return _m;
		}
		size_t coldim () const
		{
			return _n;
		}
		size_t size () const
		{
			size_t s(0);
			for(typename Rep::const_iterator it = _matA.begin(); it != _matA.end(); ++it)
				s+= LinBox::RawVector<Element>::size(*it);
			return s;
		}

		std::istream &read (std::istream &is
				    // , const Field &F
				    , LINBOX_enum(Tag::FileFormat) format /*  = Tag::FileFormat::Detect */)
		{
			return SparseMatrixReadWriteHelper<Field, Row>::read (*this, is,format);
		}

		/// Read from matrix market format
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


		std::ostream &write (std::ostream &os
				     // , const Field &F
				     , LINBOX_enum(Tag::FileFormat) format /* = Tag::FileFormat::Pretty */) const
		{
			return SparseMatrixWriteHelper<Field, Row, myTrait>::write (*this, os, format);
		}

		/// Write in matrix market format
		std::ostream &write (std::ostream &os) const
		{
			// typedef SparseMatrixBase<Element, _Row> SMB;
			writeMMCoordHeader(os, *this, this->size(), "SparseMatrix");
			return this->write(os, Tag::FileFormat::OneBased);
		}


		void           setEntry (size_t i, size_t j, const Element &value);
		Element       &refEntry (size_t i, size_t j);
		const Element &getEntry (size_t i, size_t j) const;
		Element       &getEntry (Element &x, size_t i, size_t j) const
		{
			return x = getEntry (i, j);
		}

		typedef typename Rep::iterator RowIterator;
		typedef typename Rep::const_iterator ConstRowIterator;

		ConstRowIterator rowBegin () const
		{
			return _matA.begin ();
		}
		ConstRowIterator rowEnd () const
		{
			return _matA.end ();
		}
		RowIterator rowBegin ()
		{
			return _matA.begin ();
		}
		RowIterator rowEnd ()
		{
			return _matA.end ();
		}

		template <class RepIterator, class RowEltIterator, class _I_Element>
		class _Iterator {
		public:
			typedef _I_Element value_type;

			_Iterator (const RepIterator &i, const RowEltIterator &j, const RepIterator &A_end) :
				_i (i), _j (j), _A_end (A_end)
			{
				if( _i == _A_end ) return;
				while( _j == _i->second.end() ) {
					if( ++_i == _A_end ) return;
					_j = _i->second.begin();
				}
			}

			_Iterator (const _Iterator &iter) :
				_i (iter._i), _j (iter._j), _A_end (iter._A_end)
			{}

			_Iterator () {}

			_Iterator &operator = (const _Iterator &iter)
			{
				_i = iter._i;
				_j = iter._j;
				_A_end = iter._A_end;

				return *this;
			}

			bool operator == (const _Iterator &i) const
			{
				return (_i == i._i) && (_j == i._j);
			}

			bool operator != (const _Iterator &i) const
			{
				return (_i != i._i) || (_j != i._j);
			}

			_Iterator &operator ++ ()
			{
				++_j;
				while( _j == _i->second.end() ) {
					if( ++_i == _A_end ) return *this;
					_j = _i->second.begin();
				}
				return *this;
			}

			_Iterator operator ++ (int)
			{
				_Iterator tmp = *this;
				++(*this);
				return tmp;
			}

			_Iterator &operator -- ()
			{
				while (_j == _i->second.begin ())
					_j = (--_i)->second.end ();
				--_j;
				return *this;
			}

			_Iterator operator -- (int)
			{
				_Iterator tmp = *this;
				--(*this);
				return tmp;
			}

			value_type &operator * ()
			{
				return static_cast<value_type&>(*_j);
			}
			value_type *operator -> ()
			{
				return &(*_j);
			}

		private:
			RepIterator _i;
			RowEltIterator _j;
			RepIterator _A_end;
		};

		typedef _Iterator<typename Rep::iterator, typename Row::second_type::iterator, Element> Iterator;
		typedef _Iterator<typename Rep::const_iterator, typename Row::second_type::const_iterator, const Element> ConstIterator;

		Iterator Begin ()
		{
			return Iterator (_matA.begin (), _matA.front ().second.begin (), _matA.end ());
		}
		Iterator End ()
		{
			return Iterator (_matA.end (), _matA.back ().second.end (), _matA.end ());
		}
		ConstIterator Begin () const
		{
			return ConstIterator (_matA.begin (), _matA.front ().second.begin (), _matA.end ());
		}
		ConstIterator End () const
		{
			return ConstIterator (_matA.end (), _matA.back ().second.end (), _matA.end ());
		}

		/* Generic trait for iterators without type */
		template<typename U>
		struct IteratorValueType {
			typedef typename U::value_type value_type;
		};

		template<typename X>
		struct IteratorValueType<const X*> {
			typedef X value_type;
		};


		template <class RepIterator, class RowIdxIterator>
		class _IndexedIterator {
		public:
			typedef typename IteratorValueType<RepIterator>::value_type PairValue ;
			typedef typename PairValue::second_type::value_type value_type;

			typedef _IndexedIterator<RepIterator,RowIdxIterator> Self_t;

			// typedef typename IteratorValueType< RepIterator >::second_type::value_type value_type;

			// Dan Roche 7-6-05 Fixed a seg fault this code was causing
			_IndexedIterator (size_t idx, const RepIterator &i,
					  const RowIdxIterator &j, const RepIterator &A_end) :
				_i (i), _j (j), _A_end (A_end), _r_index (idx), _c_index(0), _value_index(0)
			{
				if( _i == _A_end ) return;
				while( _j == _i->first.end() ) {
					if( ++_i == _A_end ) return;
					++_r_index;
					_j = _i->first.begin();
				}

				_c_index = *_j;
				_value_index = (size_t)( _j-_i->first.begin());
			}

			_IndexedIterator (const _IndexedIterator &iter) :
				_i (iter._i), _j (iter._j),
				_A_end (iter._A_end), _r_index (iter._r_index),
				_c_index (iter._c_index), _value_index( iter._value_index )
			{}

			_IndexedIterator ()
			{}

			_IndexedIterator &operator = (const _IndexedIterator &iter)
			{
				_A_end = iter._A_end;
				_i = iter._i;
				_j = iter._j;
				_r_index = iter._r_index;
				_c_index = iter._c_index;
				_value_index = iter._value_index;

				return *this;
			}

			bool operator == (const _IndexedIterator &i) const
			{
				return (_i == i._i) && (_j == i._j);
			}

			bool operator != (const _IndexedIterator &i) const
			{
				return (_i != i._i) || (_j != i._j);
			}

			_IndexedIterator &operator ++ ()
			{
				if(_j != _i->first.end ()) {
					++_j ;
					++_value_index;
				}
				while(_j == _i->first.end ()) {
					++_r_index;
					if (++_i == _A_end) return *this;
					_j = _i->first.begin ();
					_value_index = 0;
				}
				_c_index = *_j;

				return *this;
			}

			_IndexedIterator operator ++ (int)
			{
				_IndexedIterator tmp = *this;
				++(*this);
				return tmp;
			}

			_IndexedIterator &operator -- ()
			{
				while (_j == _i->first.begin ()) {
					_j = (--_i)->first.end ();
					_value_index = _i->first.size();
					--_r_index;
				}

				--_j;
				--_value_index;
				_c_index = *_j;
				return *this;
			}

			_IndexedIterator operator -- (int)
			{
				_IndexedIterator tmp = *this;
				--(*this);
				return tmp;
			}


			// JGD 26.11.2012
			// Since siome compliers would not choose it even though they are
			// called via a ConstIterator, const version is removed,
			// call to const is now only explicit
			// via call to "value()" below instead
			//         const value_type &operator * () const
			// 		{
			// 			return *(_i->second.begin () + _c_index);
			// 		}

			value_type &operator * ()
			{
				return (_i->second)[_value_index];
			}
#if 0
			value_type *operator -> ()
			{
				return &(*(_i->second.begin () + _c_index));
			}
			const value_type *operator -> () const
			{
				return &(*(_i->second.begin () + _c_index));
			}
#endif

			size_t rowIndex () const
			{
				return _r_index;
			}
			size_t colIndex () const
			{
				return _c_index;
			}
			const value_type &value () const
			{
				return (const value_type&)(_i->second)[_value_index];
			}

		private:
			RepIterator _i;
			RowIdxIterator _j;
			RepIterator _A_end;

			mutable size_t _r_index;
			mutable size_t _c_index;
			mutable size_t _value_index;
		};

		typedef _IndexedIterator<typename Rep::iterator, typename Row::first_type::iterator> IndexedIterator;
		typedef _IndexedIterator<typename Rep::const_iterator, typename Row::first_type::const_iterator> ConstIndexedIterator;

		IndexedIterator IndexedBegin ()
		{
			return IndexedIterator (0, _matA.begin (), _matA.front ().first.begin (), _matA.end ());
		}
		IndexedIterator IndexedEnd ()
		{
			return IndexedIterator (_m, _matA.end (), _matA.back ().first.end (), _matA.end ());
		}
		ConstIndexedIterator IndexedBegin () const
		{
			return ConstIndexedIterator (0, _matA.begin (), _matA.front ().first.begin (), _matA.end ());
		}
		ConstIndexedIterator IndexedEnd () const
		{
			return ConstIndexedIterator (_m, _matA.end (), _matA.back ().first.end (), _matA.end ());
		}

		Row &getRow (size_t i) {
			return _matA[i];
		}
		Row &operator [] (size_t i) {
			return _matA[i];
		}
		ConstRow &operator [] (size_t i) const
		{
			return _matA[i];
		}

		template <class Vector> Vector &columnDensity (Vector &v) const;
		SparseMatrix &transpose (SparseMatrix &AT) const;

		const Field & field() const
		{
			return _field ;
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

		const Field & _field;
		Rep               _matA;
		size_t            _m;
		size_t            _n;

		MatrixDomain<Field>       _MD; // Matrix domain for matrix operations
		TransposeMatrix<SparseMatrix<_Field, _Row> > _AT;

		// template<class F, class R, class T> friend class SparseMatrix;
	};

} // namespace LinBox

#include "linbox/matrix/sparse-parallel-vector.inl"

#endif // __LINBOX_matrix_sparse_parallel_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s