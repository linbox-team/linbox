/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/matrix/blas-matrix.h
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *
 * Written by :
 *               Pascal Giorgi  <pascal.giorgi@ens-lyon.fr>
 *               Clément Pernet <clement.pernet@imag.fr>
 *               Brice Boyer    <bboyer@imag.fr>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/*!@internal
 * @file matrix/blas-matrix.inl
 * @ingroup matrix
 * A \c BlasMatrix<\c _element > represents a matrix as an array of
 * <code>_element</code>s.
 */

#ifndef __LINBOX_blas_submatrix_INL
#define __LINBOX_blas_submatrix_INL

/////////////////
//   PRIVATE   //
/////////////////

namespace LinBox
{
}

//////////////////
// CONSTRUCTORS //
//////////////////

namespace LinBox
{
	template <class _Element>
	BlasSubmatrix<_Element>::BlasSubmatrix () :
			_M(NULL),_row(0),_col(0),_r0(0),_c0(0),_stride(0),_off(0)
		{
#ifndef NDEBU
			std::cerr << "***Warning*** you are creating a NULL BlasSubmatrix. You are stuck with it now." << std::endl;
#endif
		}

	template <class _Element>
	BlasSubmatrix<_Element>::BlasSubmatrix (const BlasMatrix<_Element> &Mat,
						size_t row,
						size_t col,
						size_t Rowdim,
						size_t Coldim) :
		_M (const_cast<BlasMatrix<_Element>*>(&Mat)),
		_row (Rowdim), _col(Coldim),
		_r0(row),_c0(col),
		_stride(Mat.coldim()),_off(row*_stride+col)
	{
		linbox_check ( _r0  <= Mat.rowdim() ); // allow for NULL matrix
		linbox_check ( _c0  <= Mat.coldim() );
		linbox_check ( _row <= Mat.rowdim() );
		linbox_check ( _off <= Mat.rowdim()*Mat.coldim() );
		linbox_check ( _col <= Mat.coldim() );
	}


	template <class _Element>
	BlasSubmatrix<_Element>::BlasSubmatrix (const BlasMatrix<_Element> &Mat) :
		_M(&Mat),
		_row(Mat.rowdim()), _col(Mat.coldim()),
		_r0(0), _c0(0),
		_stride(Mat.coldim()),_off(0)
	{}



	template <class _Element>
	BlasSubmatrix<_Element>::BlasSubmatrix (const BlasSubmatrix<_Element> &SM,
						size_t row,
						size_t col,
						size_t Rowdim,
						size_t Coldim) :
		_M (SM._M),
		_row ( Rowdim ),      _col ( Coldim ) ,
		_r0  ( SM._r0 + row ), _c0 ( SM._c0 + col ),
		_stride(SM._stride),_off(SM._off+(row*_stride+col))
	{
		linbox_check ( _r0  <= SM._M->rowdim() ); // allow for NULL matrix
		linbox_check ( _c0  <= SM._stride );
		linbox_check ( _row <= SM._M->rowdim() );
		linbox_check ( _col <= SM._stride );
		linbox_check ( _off <= SM._M->rowdim()*(SM._M->coldim()+1) );
	}

	template <class _Element>
	BlasSubmatrix<_Element>::BlasSubmatrix (const BlasSubmatrix<_Element> &SM) :
		_M (SM._M),
		_row ( SM._row ),  _col ( SM._col ) ,
		_r0  ( SM._r0 ),    _c0 (SM._c0 ),
		_stride(SM._stride),_off(SM._off)
	{
		linbox_check ( _r0  <=  SM._M->rowdim() ); // allow for NULL matrix
		linbox_check ( _c0  <=  SM._stride );
		linbox_check ( _row <= SM._M->rowdim() );
		linbox_check ( _col <= SM._stride );
		linbox_check ( _off <= SM._M->rowdim()*(SM._M->coldim()+1) );
	}

	template <class _Element>
	BlasSubmatrix<_Element>& BlasSubmatrix<_Element>::operator=(const BlasSubmatrix<_Element> &SM)
	{
		_M   = SM._M  ;
		_r0  = SM._r0 ;
		_row = SM._row;
		_c0  = SM._c0 ;
		_col = SM._col;
		_stride = SM._stride;
		_off = SM._off ;

		return *this;
	}

} // LinBox

///////////////////
//   ITERATORS   //
///////////////////

namespace LinBox
{


	/*! Raw Iterators.
	 * @ingroup iterators
	 *
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */
	template <class _Element>
	class BlasSubmatrix<_Element>::Iterator {
	public:
		Iterator (){}

		/*! @internal
		 * @brief NO DOC
		 */
		Iterator (const typename BlasMatrix<_Element>::Iterator& cur,
			     const size_t c_dim,
			     const size_t stride,
			     const size_t c_idx) :
			_cur (cur), _c_dim (c_dim), _stride(stride), _c_idx (c_idx)
		{}

		/*! @internal
		 * @brief copy operator.
		 * @param r Iterator to copy.
		 */
		Iterator& operator = (const Iterator& r)
		{
			_cur    = r._cur;
			_c_dim  = r._c_dim;
			_stride = r._stride;
			_c_idx  = r._c_idx;
			return *this;
		}

		/*! @internal
		 * increment.
		 * ??
		 */
		Iterator& operator ++()
		{
			if (_c_idx < _c_dim - 1){
				++_cur; ++_c_idx;
			}
			else {
				_cur = _cur + _stride - _c_dim + 1;
				_c_idx = 0;
			}

			return *this;
		}

		/*! @internal
		 * increment.
		 * ??
		 */
		Iterator& operator++ (int)
		{
			return this->operator++ ();
		}


		/*! @internal
		 * @brief  operator !=.
		 * @param r Iterator to test inequaltity from.
		 */
		bool operator != (const Iterator& r) const
		{
			return (_cur != r._cur || _c_dim != r._c_dim) || (_stride != r._stride) || (_c_idx != r._c_idx);
		}

		//! @internal operator *.
		_Element& operator * ()
		{
			return *_cur;
		}

		//! @internal operator *.
		const _Element& operator * () const
		{
			return *_cur;
		}

	protected:
		typename BlasMatrix<_Element>::Iterator _cur;
		size_t _c_dim;
		size_t _stride;
		size_t _c_idx;
	};

	/*! Raw Iterators (const version).
	 * @ingroup iterators
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */
	template <class _Element>
	class BlasSubmatrix<_Element>::ConstIterator {
	public:
		//! @internal Null constructor
		ConstIterator (){}


		/*! @internal
		 * @brief NO DOC
		 */
		ConstIterator (const typename BlasMatrix<_Element>::ConstIterator& cur,
				  const size_t c_dim,
				  const size_t stride,
				  const size_t c_idx) :
			_cur (cur), _c_dim (c_dim), _stride(stride), _c_idx (c_idx)
		{}

		/*! @internal
		 * @brief copy operator.
		 * @param r Iterator to copy.
		 */
		ConstIterator& operator = (const ConstIterator& r)
		{
			_cur = r._cur;
			_c_dim = r._c_dim;
			_stride = r._stride;
			_c_idx = r._c_idx;
			return *this;
		}

		/*! @internal
		 * increment.
		 * ??
		 */
		ConstIterator& operator ++()
		{
			if (_c_idx < _c_dim - 1){
				++_cur; ++_c_idx;
			}
			else {
				_cur = _cur + _stride - _c_dim + 1;
				_c_idx = 0;
			}

			return *this;
		}

		/*! @internal
		 * increment.
		 * ??
		 */
		ConstIterator& operator++ (int)
		{
			return this->operator++ ();
		}

		/*! @internal
		 * @brief  operator !=.
		 * @param r Iterator to test inequaltity from.
		 */
		bool operator != (const ConstIterator& r) const
		{
			return (_cur != r._cur) || (_c_dim != r._c_dim) || (_stride != r._stride) || (_c_idx != r._c_idx);
		}

		//! @internal operator *.
		const _Element& operator * () const
		{
			return *_cur;
		}

	protected:
		typename BlasMatrix<_Element>::ConstIterator _cur;
		size_t _c_dim;
		size_t _stride;
		size_t _c_idx;
	};

#if 0
	template <class Element>
	class BlasSubmatrix<Element>::ConstIterator {
	public:
		ConstIterator (){}

		ConstIterator ( const typename BlasMatrix<Element>::ConstIterator& cur,
				   size_t cont_len,
				   size_t gap_len) :
			_beg (beg), _cur (cur), _cont_len (cont_len), _gap_len (gap_len)
		{}

		ConstIterator& operator = (const Iterator& r)
		{
			_cur = r._cur;
			_beg = r._beg;
			_cont_len = r._cont_len;
			_gap_len = r._gap_len;
			return *this;
		}

		ConstIterator& operator = (const ConstIterator& r)
		{
			_cur = r._cur;
			_beg = r._beg;
			_cont_len = r._cont_len;
			_gap_len = r._gap_len;
			return *this;
		}

		ConstIterator& operator++()
		{
			if (((_cur - _beg + 1) % _cont_len) != 0)
				++_cur;
			else
			{
				_cur = _cur + _gap_len + 1;
				_beg = _beg + _gap_len + _cont_len;
			}
			return *this;
		}

		ConstIterator operator++(int)
		{
			ConstIterator tmp = *this;
			this->operator++();
			return tmp;
		}

		bool operator != (const ConstIterator& r) const
		{
			return (_cur != r._cur) || (_beg != r._beg) || (_cont_len != r._cont_len) || (_gap_len != r._gap_len);
		}

		const Element& operator*()
		{ return *_cur; }

		Element& operator*()
		{ return *_cur; }

		const Element& operator*() const
		{ return *_cur; }

	protected:
		typename BlasMatrix<Element>::ConstIterator _beg;
		typename BlasMatrix<Element>::ConstIterator _cur;
		size_t _cont_len;
		size_t _gap_len;
	};
#endif

	template <class _Element>
	typename BlasSubmatrix<_Element>::Iterator BlasSubmatrix<_Element>::Begin ()
	{
		return Iterator (_M->Begin () + ( _off ),
				    _col, _stride, 0);
	}

	template <class _Element>
	typename BlasSubmatrix<_Element>::Iterator BlasSubmatrix<_Element>::End ()
	{
		return Iterator (_M->Begin () + ( (_row) * _stride + _off ),
				    _col, _stride, 0);
	}

	template <class _Element>
	typename BlasSubmatrix<_Element>::ConstIterator BlasSubmatrix<_Element>::Begin () const
	{
		return ConstIterator (_M->Begin () + ( _off ),
					 _col, _stride, 0);
	}

	template <class _Element>
	typename BlasSubmatrix<_Element>::ConstIterator BlasSubmatrix<_Element>::End () const
	{
		return ConstIterator (_M->Begin () + ( (_row) * _stride + _off ),
					 _col, _stride, 0);
	}

#if 0
	template <class Element>
	typename BlasSubmatrix<Element>::ConstIterator BlasSubmatrix<Element>::Begin () const
	{
		return ConstIterator (_M->Begin () + ( _off ),
					 _M->Begin () + ( _off ),
					 _col, _stride - _col);
	}

	template <class Element>
	typename BlasSubmatrix<Element>::ConstIterator BlasSubmatrix<Element>::End () const
	{
		return ConstIterator (_M->Begin () + ( (_row) * _stride + _off ),
					 _M->Begin () + ( (_row) * _stride + _off ),
					 _col, _stride - _col);
	}
#endif

	/*! Raw Indexed Iterator.
	 * @ingroup iterators
	 *
	 * Like the raw iterator, the indexed iterator is a method for
	 * accessing all entries in the matrix in some unspecified order.
	 * At each position of the the indexed iterator, it also provides
	 * the row and column indices of the currently referenced entry.
	 * This is provided through it's \c rowIndex() and \c colIndex() functions.
	 */
	template <class _Element>
	class BlasSubmatrix<_Element>::IndexedIterator {
	public:
		IndexedIterator (){}

		IndexedIterator (const typename BlasMatrix<_Element>::Iterator& cur,
				    size_t c_dim,
				    size_t stride,
				    size_t r_idx,
				    size_t c_idx) :
			_cur (cur), _c_dim (c_dim), _stride (stride), _r_idx(r_idx), _c_idx (c_idx)
		{}

		IndexedIterator& operator = (const IndexedIterator& r)
		{
			_cur = r._cur;
			_stride = r._stride;
			_c_dim = r._c_dim;
			_r_idx = r._r_idx;
			_c_idx = r._c_idx;
			return *this;
		}

		IndexedIterator& operator++()
		{
			if (_c_idx < _c_dim - 1){
				++_c_idx;
				++_cur;
			}
			else
			{
				_cur = _cur + _stride - _c_dim + 1;
				_c_idx = 0;
				++_r_idx;
			}
			return *this;
		}

		IndexedIterator& operator--()
		{
			if (_c_idx > 0){
				--_c_idx;
				--_cur;
			}
			else
			{
				_cur = _cur - _stride + _c_dim -1;
				_c_idx = 0;
				--_r_idx;
			}
			return *this;
		}

		IndexedIterator operator++(int)
		{
			IndexedIterator tmp = *this;
			this->operator++();
			return tmp;
		}

		IndexedIterator operator--(int)
		{
			IndexedIterator tmp = *this;
			this->operator--();
			return tmp;
		}

		bool operator != (const IndexedIterator& r) const
		{
			return ((_c_idx != r._c_idx) || (_r_idx != r._r_idx) ||(_stride != r._stride) || (_c_dim != r._c_dim) );
		}

		const _Element& operator*() const {return *_cur;}

		_Element& operator*() {return *_cur;}

		size_t rowIndex () const { return _r_idx; }

		size_t colIndex () const { return _c_idx; }

		const _Element& value () const {return *_cur;}

	protected:
		typename BlasMatrix<_Element>::Iterator _cur;
		size_t _stride;
		size_t _c_dim;
		size_t _r_idx;
		size_t _c_idx;
	};

	template <class _Element>
	typename BlasSubmatrix<_Element>::IndexedIterator BlasSubmatrix<_Element>::IndexedBegin ()
	{
		return IndexedIterator (_M->Begin () + ( (_off) ),
					   _col , _stride, 0, 0);
	}

	template <class _Element>
	typename BlasSubmatrix<_Element>::IndexedIterator BlasSubmatrix<_Element>::IndexedEnd ()
	{
		return IndexedIterator (_M->Begin () + ( (_row) * _stride + (_col+_off) ),
					   _col, _stride, _row-1, _col-1);
	}

	/*! Raw Indexed Iterator (const version).
	 * @ingroup iterators
	 *
	 * Like the raw iterator, the indexed iterator is a method for
	 * accessing all entries in the matrix in some unspecified order.
	 * At each position of the the indexed iterator, it also provides
	 * the row and column indices of the currently referenced entry.
	 * This is provided through it's \c rowIndex() and \c colIndex() functions.
	 */
	template <class _Element>
	class BlasSubmatrix<_Element>::ConstIndexedIterator {
	public:
		ConstIndexedIterator (){}

		ConstIndexedIterator (const typename BlasMatrix<_Element>::ConstIterator& cur,
					 size_t c_dim,
					 size_t stride,
					 size_t r_idx,
					 size_t c_idx) :
			_cur (cur), _stride (stride), _c_dim (c_dim), _r_idx(r_idx), _c_idx (c_idx)
		{}

		ConstIndexedIterator& operator = (const IndexedIterator& r)
		{
			_cur = r._cur;
			_stride = r._stride;
			_c_dim = r._c_dim;
			_r_idx = r._r_idx;
			_c_idx = r._c_idx;
			return *this;
		}

		ConstIndexedIterator& operator = (const ConstIndexedIterator& r)
		{
			_cur = r._cur;
			_stride = r._stride;
			_c_dim = r._c_dim;
			_r_idx = r._r_idx;
			_c_idx = r._c_idx;
			return *this;
		}

		ConstIndexedIterator& operator++()
		{
			if (_c_idx < _c_dim - 1){
				++_c_idx;
				++_cur;
			}
			else
			{
				_cur = _cur + _stride - _c_dim +1;
				_c_idx = 0;
				++_r_idx;
			}
			return *this;
		}

		IndexedIterator& operator--()
		{
			if (_c_idx > 0){
				--_c_idx;
				--_cur;
			}
			else
			{
				_cur = _cur - _stride + _c_dim -1;
				_c_idx = 0;
				--_r_idx;
			}
			return *this;
		}

		ConstIndexedIterator operator++(int)
		{
			ConstIndexedIterator tmp = *this;
			this->operator++();
			return tmp;
		}

		ConstIndexedIterator operator--(int)
		{
			ConstIndexedIterator tmp = *this;
			this->operator--();
			return tmp;
		}

		size_t rowIndex () const
		{
			return _r_idx;
		}

		size_t colIndex () const
		{
			return _c_idx;
		}

		bool operator != (const ConstIndexedIterator& r) const
		{
			return ((_c_idx != r._c_idx) || (_r_idx != r._r_idx) ||(_stride != r._stride) || (_c_dim != r._c_dim) );
		}

		const _Element& operator*() const
		{
			return *_cur;
		}


		friend std::ostream& operator<<(std::ostream& out, const ConstIndexedIterator m)
		{
			return out /* << m._cur << ' ' */
			<< m._stride << ' '
			<< m._c_dim << ' '
			<< m._r_idx << ' '
			<< m._c_idx;
		}

		const Element & value() const
		{
			return this->operator*();

		}

	protected:
		typename BlasMatrix<_Element>::ConstIterator _cur;
		size_t _stride;
		size_t _c_dim;
		size_t _r_idx;
		size_t _c_idx;
	};

	template <class _Element>
	typename BlasSubmatrix<_Element>::ConstIndexedIterator BlasSubmatrix<_Element>::IndexedBegin () const
	{
		return ConstIndexedIterator (_M->Begin () + ( _off ),
						_row, _stride, 0, 0);
	}

	template <class _Element>
	typename BlasSubmatrix<_Element>::ConstIndexedIterator BlasSubmatrix<_Element>::IndexedEnd () const
	{
		return ConstIndexedIterator (_M->Begin () + ( (_row) * _stride + (_off+_col) ),
						_col, _stride, _row-1, _col-1);
	}

	////////
	template <class _Element>
	typename BlasSubmatrix<_Element>::RowIterator BlasSubmatrix<_Element>::rowBegin ()
	{
		return RowIterator (_M->Begin () + ( _off  ),
				    _col, _stride);
	}

	template <class _Element>
	typename BlasSubmatrix<_Element>::RowIterator BlasSubmatrix<_Element>::rowEnd ()
	{
		return RowIterator (_M->Begin () + ( (_row) * _stride + _off ),
				    _col, _stride);
	}

	template <class _Element>
	typename BlasSubmatrix<_Element>::ConstRowIterator BlasSubmatrix<_Element>::rowBegin () const
	{
		return ConstRowIterator (_M->Begin () + ( _off ),
					 _col, _stride);
	}

	template <class _Element>
	typename BlasSubmatrix<_Element>::ConstRowIterator BlasSubmatrix<_Element>::rowEnd () const
	{
		return ConstRowIterator (_M->Begin () + ( (_row) * _stride + _off ),
					 _col, _stride);
	}

	template <class _Element>
	typename BlasSubmatrix<_Element>::ColIterator BlasSubmatrix<_Element>::colBegin ()
	{
		return ColIterator (_M->Begin () + ( _off ),
				    _stride, _row);
	}

	template <class _Element>
	typename BlasSubmatrix<_Element>::ColIterator BlasSubmatrix<_Element>::colEnd ()
	{
		return ColIterator (_M->Begin () + ( (_col) + _off ),
				    _stride, _row);
	}

	template <class _Element>
	typename BlasSubmatrix<_Element>::ConstColIterator BlasSubmatrix<_Element>::colBegin () const
	{
		return ConstColIterator (_M->Begin () + ( _off ),
					 _stride, _row);
	}

	template <class _Element>
	typename BlasSubmatrix<_Element>::ConstColIterator BlasSubmatrix<_Element>::colEnd () const
	{
		return ConstColIterator (_M->Begin () + ( (_col) + _off ),
					 _stride, _row);
	}

	/*  operators */
	template <class Element>
	typename BlasSubmatrix<Element>::Row BlasSubmatrix<Element>::operator[] (size_t i)
	{
		return Row (_M.Begin () + (_r0+i) * _stride, _M.Begin () + ((_r0+i) * _stride + _stride) );
	}

	template <class Element>
	typename BlasSubmatrix<Element>::ConstRow BlasSubmatrix<Element>::operator[] (size_t i) const
	{
		return Row (_M.Begin () + (_r0+i) * _stride, _M.Begin () + ((_r0+i) * _stride + _stride) );
	}

}

///////////////////
//      I/O      //
///////////////////

namespace LinBox
{
	template <class _Element>
	template <class Field>
	std::istream& BlasSubmatrix<_Element>::read (std::istream &file, const Field& field)
	{
		Iterator p;

		for (p = Begin (); p != End (); ++p) {
			// each entry is seperated by one space.
			file.ignore (1);
			field.read (file, *p);
		}

		return file;
	}

	template <class _Element>
	template <class Field>
	std::ostream &BlasSubmatrix<_Element>::write (std::ostream &os, const Field& field,
						       bool mapleFormat) const
	{
		ConstRowIterator p;

		//!@todo Integer width ?

		typename ConstRow::const_iterator pe;


		if (mapleFormat) os << "Matrix( " << _row << ',' << _col << ",[" ;

		for (p = rowBegin (); p != rowEnd (); ++p) {
			if (mapleFormat && (p != rowBegin()))
				os << ',';
			if (mapleFormat) os << "[";

			for (pe = p->begin (); pe != p->end (); ++pe) {
				if (mapleFormat && (pe != p->begin())) os << ',';
				/*!  @warning
				* matrix base does not provide this field(), maybe should?
				* _M.field ().write (os, *pe);
				* os << *pe;
				* fixed by using extra field
				*/

				field.write (os, *pe);
				os << " ";
			}

			if (!mapleFormat)
				os << std::endl;
			else os << ']';
		}

		if (mapleFormat) os << "])";
		return os;
	}

	template <class _Element>
	std::ostream &BlasSubmatrix<_Element>::write (std::ostream &os, bool mapleFormat) const
	{
		ConstRowIterator p;



		typename ConstRow::const_iterator pe;

		if (mapleFormat) os << "Matrix( " << _row << ',' << _col << ",[" ;

		for (p = rowBegin (); p != rowEnd (); ++p) {
			if (mapleFormat && (p != rowBegin()))
				os << ',';
			if (mapleFormat) os << "[";

			for (pe = p->begin (); pe != p->end (); ++pe) {
				if (mapleFormat && (pe != p->begin())) os << ',';

				os << *pe;
				os << " ";
			}

			if (!mapleFormat)
				os << std::endl;
			else os << ']';
		}

		if (mapleFormat) os << "])";
		return os;
	}


} // LinBox

//////////////////
//  DIMENSIONS  //
//////////////////
namespace LinBox
{
	template <class Element>
	size_t BlasSubmatrix<Element>::rowdim() const
	{
		return _row ;
	}

	template <class Element>
	size_t BlasSubmatrix<Element>::coldim() const
	{
		return _col ;
	}

	template <class Element>
	size_t  BlasSubmatrix<Element>::getStride() const
	{
		return _stride;
	}


}

//////////////////
//   ELEMENTS   //
//////////////////


namespace LinBox
{
	template <class Element>
	void BlasSubmatrix<Element>::setEntry (size_t i, size_t j, const Element &a_ij)
	{
		_M->setEntry (_r0 + i, _c0 + j, a_ij);
	}

	template <class Element>
	Element &BlasSubmatrix<Element>::refEntry (size_t i, size_t j)
	{
		return _M->refEntry ( _r0+i, _c0+j );
	}

	template <class Element>
	const Element &BlasSubmatrix<Element>::getEntry (size_t i, size_t j) const
	{
		return _M->getEntry ( _r0+i , _c0+j );
	}

	template <class Element>
	Element &BlasSubmatrix<Element>::getEntry (Element &x, size_t i, size_t j) const
	{
		return _M->getEntry (x, _r0+i , _c0+j );
	}

}
#endif // __LINBOX_blas_submatrix_INL
