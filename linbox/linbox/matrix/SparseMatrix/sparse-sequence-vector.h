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
 *   - Made SparseMatrixGeneric parameterized only on the element type
 *   - New read/write implementations for SparseMatrixGeneric, supporting multiple
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

/** @file linbox/matrix/SparseMatrix/sparse-sequence-vector.h
 * @brief
 */

#ifndef __LINBOX_matrix_sparsematrix_sparse_sequence_vector_H
#define __LINBOX_matrix_sparsematrix_sparse_sequence_vector_H


namespace LinBox { namespace Protected {


	/* Specialization for sparse sequence vectors */
	template <class _Field, class _Row>
	class SparseMatrixGeneric<_Field, _Row, VectorCategories::SparseSequenceVectorTag > {
	public:

		typedef         _Field          Field;
		typedef typename Field::Element Element;

		typedef      _Row           Row;
		typedef const Row           ConstRow;
		typedef _SP_BB_VECTOR_<Row> Rep;
		typedef VectorCategories::SparseSequenceVectorTag myTrait;
		typedef SparseMatrixGeneric<_Field, _Row, myTrait> Self_t;

#ifdef __LINBOX_PARALLEL
		BB_list_list sub_list;
#endif


		template<typename _Tp1, typename _R1 = typename Rebind<_Row,_Tp1>::other >
		struct rebind {
			typedef SparseMatrixGeneric<_Tp1, _R1> other;

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
		SparseMatrixGeneric (const SparseMatrixGeneric<_Tp1, _Rw1, myTrait> &Mat, const Field& F) :
			_field (F),
			_MD (F), _AT (*this)
			, _matA (Mat.rowdim()), _m (Mat.rowdim()), _n (Mat.coldim())
		{
			typename SparseMatrixGeneric<_Tp1,_Rw1, myTrait>::template rebind<Field,_Row>()(*this, Mat);
		}

		SparseMatrixGeneric (const Field & F, size_t m, size_t n) :
			_field(F),
			_MD(F),_AT(*this),
			_matA (m), _m (m), _n (n)
		{}

		SparseMatrixGeneric (const Field & F) :
			_field (F),
			_MD(F),_AT(*this),
			_matA(0), _m(0), _n(0)
		{};


		SparseMatrixGeneric (const SparseMatrixGeneric<Field, Row> &A) :
			_field(A.field()),
			_MD(A.field()),_AT(*this),
			_matA (A._matA), _m (A._m), _n (A._n)
		{}

		template<class VectorType>
		SparseMatrixGeneric (const SparseMatrixGeneric<Field, VectorType> &A) :
			_field(A.field()),
			_MD(A.field()),_AT(*this),
			_matA(A.rowdim()), _m (A.rowdim()), _n (A.coldim())
		{
			typename Rep::iterator meit = this->_matA.begin();
			typename SparseMatrixGeneric<Field, VectorType>::Rep::const_iterator copit = A.getRep().begin();
			for( ; meit != this->_matA.end(); ++meit, ++copit)
				LinBox::RawVector<Element>::convert(*meit, *copit);
		}

		/** Constructor from a MatrixStream
		*/
		SparseMatrixGeneric ( MatrixStream<Field>& ms );

		template<class VectStream>
		SparseMatrixGeneric (const Field &F, VectStream &stream) :
			_field (F), _MD (F), _AT (*this)
			, _matA (stream.size()), _m (stream.size()), _n (stream.dim())
		{
			typename Self_t::RowIterator i;

			for (i = Self_t::rowBegin (); i != Self_t::rowEnd (); ++i)
				stream >> *i;
		}

		~SparseMatrixGeneric () {}

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
				    , LINBOX_enum(Tag::FileFormat) format  = Tag::FileFormat::Detect )
		{
			return SparseMatrixReadHelper<Self_t>::read (*this, is
									      , format);
		}

		std::ostream &write (std::ostream &os
				     , LINBOX_enum(Tag::FileFormat) format /* = Tag::FileFormat::Pretty*/) const
		{
			return SparseMatrixWriteHelper<Self_t>::write (*this, os, format);
		}

		/// Write in matrix market format
		std::ostream &write (std::ostream &os) const
		{
			writeMMCoordHeader(os, *this, this->size(), "SparseMatrixGeneric");
			return this->write(os, Tag::FileFormat::OneBased);
		}

		void finalize(){}

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

		template <class RepIterator, class RowIterator, class _I_Element>
		class _Iterator {
		public:
			typedef _I_Element value_type;

			_Iterator (const RepIterator &i, const RowIterator &j, const RepIterator &A_end) :
				_i (i), _j (j), _A_end (A_end)
			{
				if( _i == _A_end ) return;
				while ( _j == _i->end () ) {
					if (++_i == _A_end) return;
					_j = _i->begin ();
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
				while( _j == _i->end ()) {
					if (++_i == _A_end) return *this;
					_j = _i->begin ();
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
				while (_j == _i->begin ())
					_j = (--_i)->end ();
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
				return _j->second;
			}
			// Dan Roche 2005-7-7 I believe this was a memory leak.
			value_type *operator -> ()
			{
				return &(_j->second);
			}

			const value_type &operator*() const
			{
				return _j->second;
			}

			const value_type *operator-> () const
			{
				return &(_j->second);
			}

		private:
			RepIterator _i;
			RowIterator _j;
			RepIterator _A_end;
		};

		typedef _Iterator<typename Rep::iterator, typename Row::iterator, Element> Iterator;
		typedef _Iterator<typename Rep::const_iterator, typename Row::const_iterator, const Element> ConstIterator;

		Iterator Begin ()
		{
			return Iterator (_matA.begin (), _matA.front ().begin (), _matA.end ());
		}
		Iterator End ()
		{
			return Iterator (_matA.end (), _matA.back ().end (), _matA.end ());
		}
		ConstIterator Begin () const
		{
			return ConstIterator (_matA.begin (), _matA.front ().begin (), _matA.end ());
		}
		ConstIterator End () const
		{
			return ConstIterator (_matA.end (), _matA.back ().end (), _matA.end ());
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

		/* Generic trait for iterators without type */



		template <class RepIterator, class RowIdxIterator>
		class _IndexedIterator {
		public:
#if 0
			typedef std::pair<size_t, size_t> value_type;
			typedef typename RowIdxIterator/*::value_type*/::second_type value_type;
			typedef typename RowIdxIterator::value_type::second_type value_type;
#endif
			typedef typename IteratorValueType< RowIdxIterator >::value_type::second_type value_type;

			_IndexedIterator (size_t idx, const RepIterator &i, const RowIdxIterator &j, const RepIterator &A_end) :
				_i (i), _j (j), _A_end (A_end), _r_index (idx)
			{
				if( _i == _A_end ) return;
				while(_j == _i->end ()) {
					++_r_index;
					if (++_i == _A_end) return;
					_j = _i->begin ();
				}
				_c_index =_j->first;
			}

			_IndexedIterator (const _IndexedIterator &iter) :
				_i (iter._i), _j (iter._j), _A_end (iter._A_end), _r_index (iter._r_index), _c_index (iter._c_index)
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
				++_j;
				while(_j == _i->end ()){
					++_r_index;
					if (++_i == _A_end) return *this;
					_j = _i->begin ();
				}
				_c_index = _j->first;
#if 0
				if (++_j == _i->end ()) {
					if (++_i != _A_end) {
						_j = _i->begin ();
						++_r_index;
					}
				}

				_c_index = _j->first;
#endif

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
				while (_j == _i->begin ()) {
					_j = (--_i)->end ();
					--_r_index;
				}

				--_j;
				_c_index = _j->first;
				return *this;
			}

			_IndexedIterator operator -- (int)
			{
				_IndexedIterator tmp = *this;
				--(*this);
				return tmp;
			}

			value_type &operator * ()
			{
				return *const_cast<value_type*> (&(_j->second));
				// Ugh.  This hack is because ConstIndexedIterator is not right. -bds

			}
			const value_type &operator * () const
			{
				return _j->second;
			}
			value_type *operator -> ()
			{
				return &(_j->second);
			}
			const value_type *operator -> () const
			{
				return &(_j->second);
			}

			size_t rowIndex () const
			{
				return _r_index;
			}
			size_t colIndex () const
			{
				return _c_index;
			}
			const value_type &value() const
			{
				return _j->second;
			}

		private:
			RepIterator _i;
			RowIdxIterator _j;
			RepIterator _A_end;

			mutable size_t _r_index;
			mutable size_t _c_index;
		};

		typedef _IndexedIterator<typename Rep::iterator, typename Row::iterator> IndexedIterator;
		typedef _IndexedIterator<typename Rep::const_iterator, typename Row::const_iterator> ConstIndexedIterator;

		IndexedIterator IndexedBegin ()
		{
			return IndexedIterator (0, _matA.begin (), _matA.front ().begin (), _matA.end ());
		}
		IndexedIterator IndexedEnd ()
		{
			return IndexedIterator (_m, _matA.end (), _matA.back ().end (), _matA.end ());
		}
		ConstIndexedIterator IndexedBegin () const
		{
			return ConstIndexedIterator (0, _matA.begin (), _matA.front ().begin (), _matA.end ());
		}
		ConstIndexedIterator IndexedEnd () const
		{
			return ConstIndexedIterator (_m, _matA.end (), _matA.back ().end (), _matA.end ());
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
		SparseMatrixGeneric &transpose (SparseMatrixGeneric &AT) const;

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

		Rep & refRep()
		{
			return _matA;
		}

	protected:

		friend class SparseMatrixWriteHelper<Self_t>;
		friend class SparseMatrixReadHelper<Self_t>;

		const Field & _field;

		MatrixDomain<Field>       _MD; // Matrix domain for matrix operations
		TransposeMatrix<SparseMatrixGeneric<_Field, _Row> > _AT;
		Rep               _matA;
		size_t            _m;
		size_t            _n;

		// template<class F, class R, class T> friend class SparseMatrixGeneric;
	};

} // namespace LinBox
} // namespace Protected



#include "linbox/matrix/SparseMatrix/sparse-sequence-vector.inl"

namespace LinBox
{

	template <class _Field /*, class _Row */  >
	class SparseMatrix<_Field, SparseMatrixFormat::SparseSeq/* <_Row> */ > : public Protected::SparseMatrixGeneric<_Field,/*  _Row */ typename Vector<_Field>::SparseSeq ,VectorCategories::SparseSequenceVectorTag>
	{
	public:
		typedef VectorCategories::SparseSequenceVectorTag  myTrait ;
		typedef _Field                                       Field ; //!< Field
		typedef typename _Field::Element                   Element ; //!< Element
		typedef const Element                         constElement ; //!< const Element
		typedef typename Vector<_Field>::SparseSeq             Row ;
		typedef SparseMatrixFormat::SparseSeq              Storage ; //!< Matrix Storage Format
		typedef SparseMatrix<_Field,Storage>               Self_t ; //!< Self type
		typedef Protected::SparseMatrixGeneric<_Field,Row,myTrait >         Father_t ;

	public:
		template<class VectStream>
		SparseMatrix (const Field &F, VectStream &stream) :
			Father_t(F,stream)
		{}

		SparseMatrix(const Field & F, size_t m, size_t n) :
			Father_t(F, m, n)
		{}

		SparseMatrix(const Field & F) :
			Father_t(F)
		{}

		SparseMatrix ( MatrixStream<Field>& ms ) :
			Father_t(ms)
		{}


		using Father_t::rebind;
	} ; // SparseMatrix

	template <class Field>
	struct MatrixTraits< SparseMatrix<Field, SparseMatrixFormat::SparseSeq> >
	{
		typedef SparseMatrix<Field, SparseMatrixFormat::SparseSeq> MatrixType;
		typedef typename MatrixCategories::RowMatrixTag MatrixCategory;
	};

	template <class Field>
	struct MatrixTraits< const SparseMatrix<Field, SparseMatrixFormat::SparseSeq> >
	{
		typedef SparseMatrix<Field, SparseMatrixFormat::SparseSeq> MatrixType;
		typedef typename MatrixCategories::RowMatrixTag MatrixCategory;
	};

} // namespace LinBox

#endif // __LINBOX_matrix_sparsematrix_sparse_sequence_vector_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
