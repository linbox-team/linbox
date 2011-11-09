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
 * A \c BlasMatrix<\c _Field > represents a matrix as an array of
 * <code>_Field</code>s.
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
	template <class _Field>
	BlasSubmatrix<_Field>::BlasSubmatrix () :
			_M(NULL),_row(0),_col(0),_r0(0),_c0(0),_stride(0),_off(0)
		{
#ifndef NDEBUG
			std::cerr << "***Warning*** you are creating a NULL BlasSubmatrix. You are stuck with it now." << std::endl;
#endif
		}

	template <class _Field>
	BlasSubmatrix<_Field>::BlasSubmatrix (const BlasMatrix<_Field> &Mat,
						size_t row,
						size_t col,
						size_t Rowdim,
						size_t Coldim) :
		_M (const_cast<BlasMatrix<_Field>*>(&Mat)),
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


	template <class _Field>
	BlasSubmatrix<_Field>::BlasSubmatrix (const BlasMatrix<_Field> &Mat) :
		_M(&Mat),
		_row(Mat.rowdim()), _col(Mat.coldim()),
		_r0(0), _c0(0),
		_stride(Mat.coldim()),_off(0)
	{}



	template <class _Field>
	BlasSubmatrix<_Field>::BlasSubmatrix (const BlasSubmatrix<_Field> &SM,
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

	template <class _Field>
	BlasSubmatrix<_Field>::BlasSubmatrix (const BlasSubmatrix<_Field> &SM) :
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

	template <class _Field>
	BlasSubmatrix<_Field>& BlasSubmatrix<_Field>::operator=(const BlasSubmatrix<_Field> &SM)
	{
		if ( &SM == this)
			return *this ;

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
//      I/O      //
///////////////////

namespace LinBox
{

	//!@bug reading a submatrix should not be allowed !!
	template <class _Field>
	std::istream& BlasSubmatrix< _Field>::read (std::istream &file)
	{
#if 0
		Iterator p;
		int m,n;
		char c;
		file>>m>>n>>c;

		if (m*n < _row*_col)
			cerr<<"NOT ENOUGH ELEMENT TO READ\n";
		else {
			for (p = Begin (); p != End (); ++p) {
				integer tmp;
				file>>tmp;cout<<tmp<<endl;
				//file.ignore(1);
				_M->_F.read (file, *p);
			}
		}
#endif


		Iterator p;
		int m,n;
		char c;
		file>>m>>n>>c;
		// std::cout << m << 'x' << n << ':' << c << std::endl;
		_row = m; _col = n;

		_Field zero;
		_M->_F.init(zero,0UL);
		// resize(_row,_col);

		if ((c != 'M') && (c != 'm')) {
		for (p = Begin (); p != End (); ++p) {
				//file.ignore(1);
				_M->_F.read (file, *p);
			}

		}
		else { // sparse file format - needs fixing
			int i, j;
			while (true)
			{
				file >> i >> j;
				//file.ignore(1);
				//if (! file) break;
				if (i+j <= 0) break;
				// std::cout << i << ',' << j << ':' ;
				_M->_F.read (file, _M->refEntry[i-1, j-1]);
			}
		}

		return file;
	}

	template <class _Field>
	std::ostream &BlasSubmatrix< _Field>::write (std::ostream &os,
						       bool mapleFormat) const
	{

		ConstRowIterator p;

		if (!mapleFormat) {
			integer c;
			int wid;




			_M->_F.cardinality (c);

			if (c >0)
				wid = (int) ceil (log ((double) c) / M_LN10);
			else {
				integer tmp;
				size_t max=0;
				ConstIterator it = Begin();
				for (; it != End(); ++it){
					_M->_F.convert(tmp,*it);
					if (tmp.bitsize() > max)
						max= tmp.bitsize();
				}
				wid= (int) ceil ((double)max / M_LN10)+1;
			}

			for (p = rowBegin (); p != rowEnd ();++p) {
				typename ConstRow::const_iterator pe;

				os << "  [ ";

			for (pe = p->begin (); pe != p->end (); ++pe) {
					os.width (wid);
					/*!  @warning
					 * matrix base does not provide this field(), maybe should?
					 * _M.field ().write (os, *pe);
					 * os << *pe;
					 * fixed by using extra field
					 */

					_M->_F.write (os, *pe);
				os << " ";
			}

				os << "]" << std::endl;
			}
		}
		else {

			os << "Matrix( " << rowdim() << ',' << coldim() << ",[" ;
			for (p = rowBegin (); p != rowEnd (); ) {
				typename ConstRow::const_iterator pe;

				os << " [ ";

				for (pe = p->begin (); pe != p->end (); ) {
					_M->_F.write (os, *pe);
					++pe ;
					if (pe != p->end())
						os << ", ";
		}

				os << "]" ;
				++p ;
				if (p != rowEnd() )
					os << ',' << std::endl;;

			}
			os << "])" ;
		}
		return os;
	}



} // LinBox

//////////////////
//  DIMENSIONS  //
//////////////////

namespace LinBox
{

	template <class _Field>
	size_t BlasSubmatrix< _Field>::rowdim() const
	{
		return _row ;
	}

	template <class _Field>
	size_t BlasSubmatrix< _Field>::coldim() const
	{
		return _col ;
	}

	template <class _Field>
	size_t  BlasSubmatrix< _Field>::getStride() const
	{
		return _stride;
	}


} // LinBox

//////////////////
//   ELEMENTS   //
//////////////////

namespace LinBox
{


	template <class _Field>
	void BlasSubmatrix< _Field>::setEntry (size_t i, size_t j, const Element &a_ij)
	{
		_M->setEntry (_r0 + i, _c0 + j, a_ij);
	}

	template <class _Field>
 typename _Field:: Element &BlasSubmatrix< _Field>::refEntry (size_t i, size_t j)
	{
		return _M->refEntry ( _r0+i, _c0+j );
	}

	template <class _Field>
	const typename _Field:: Element &BlasSubmatrix< _Field>::getEntry (size_t i, size_t j) const
	{
		return _M->getEntry ( _r0+i , _c0+j );
	}

	template <class _Field>
typename _Field:: Element &BlasSubmatrix< _Field>::getEntry (Element &x, size_t i, size_t j) const
	{
		return _M->getEntry (x, _r0+i , _c0+j );
	}

}

///////////////////
// TRANSPOSE &AL //
///////////////////

namespace LinBox
{
}

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
	template <class _Field>
	class BlasSubmatrix<_Field>::Iterator {
	public:
		Iterator (){}

		/*! @internal
		 * @brief NO DOC
		 */
		Iterator (const typename BlasMatrix<_Field>::Iterator& cur,
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
		_Field& operator * ()
		{
			return *_cur;
		}

		//! @internal operator *.
		const _Field& operator * () const
		{
			return *_cur;
		}

	protected:
		typename BlasMatrix<_Field>::Iterator _cur;
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
	template <class _Field>
	class BlasSubmatrix<_Field>::ConstIterator {
	public:
		//! @internal Null constructor
		ConstIterator (){}


		/*! @internal
		 * @brief NO DOC
		 */
		ConstIterator (const typename BlasMatrix<_Field>::ConstIterator& cur,
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
		const _Field& operator * () const
		{
			return *_cur;
		}

	protected:
		typename BlasMatrix<_Field>::ConstIterator _cur;
		size_t _c_dim;
		size_t _stride;
		size_t _c_idx;
	};

#if 0
	template <class _Field>
	class BlasSubmatrix< _Field>::ConstIterator {
	public:
		ConstIterator (){}

		ConstIterator ( const typename BlasMatrix< _Field>::ConstIterator& cur,
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
			else {
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

		const _Field& operator*()
		{ return *_cur; }

	 _Field& operator*()
		{ return *_cur; }

		const _Field& operator*() const
		{ return *_cur; }

	protected:
		typename BlasMatrix< _Field>::ConstIterator _beg;
		typename BlasMatrix< _Field>::ConstIterator _cur;
		size_t _cont_len;
		size_t _gap_len;
	};
#endif

	template <class _Field>
	typename BlasSubmatrix<_Field>::Iterator BlasSubmatrix<_Field>::Begin ()
	{
		return Iterator (_M->Begin () + ( _off ),
				    _col, _stride, 0);
	}

	template <class _Field>
	typename BlasSubmatrix<_Field>::Iterator BlasSubmatrix<_Field>::End ()
	{
		return Iterator (_M->Begin () + ( (_row) * _stride + _off ),
				    _col, _stride, 0);
	}

	template <class _Field>
	typename BlasSubmatrix<_Field>::ConstIterator BlasSubmatrix<_Field>::Begin () const
	{
		return ConstIterator (_M->Begin () + ( _off ),
					 _col, _stride, 0);
	}

	template <class _Field>
	typename BlasSubmatrix<_Field>::ConstIterator BlasSubmatrix<_Field>::End () const
	{
		return ConstIterator (_M->Begin () + ( (_row) * _stride + _off ),
					 _col, _stride, 0);
	}

#if 0
	template <class _Field>
	typename BlasSubmatrix< _Field>::ConstIterator BlasSubmatrix< _Field>::Begin () const
	{
		return ConstIterator (_M->Begin () + ( _off ),
					 _M->Begin () + ( _off ),
					 _col, _stride - _col);
	}

	template <class _Field>
	typename BlasSubmatrix< _Field>::ConstIterator BlasSubmatrix< _Field>::End () const
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
	template <class _Field>
	class BlasSubmatrix<_Field>::IndexedIterator {
	public:
		IndexedIterator (){}

		IndexedIterator (const typename BlasMatrix<_Field>::Iterator& cur,
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
			else {
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
			else {
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

		const _Field& operator*() const {return *_cur;}

		_Field& operator*() {return *_cur;}

		size_t rowIndex () const { return _r_idx; }

		size_t colIndex () const { return _c_idx; }

		const _Field& value () const {return *_cur;}

	protected:
		typename BlasMatrix<_Field>::Iterator _cur;
		size_t _stride;
		size_t _c_dim;
		size_t _r_idx;
		size_t _c_idx;
	};

	template <class _Field>
	typename BlasSubmatrix<_Field>::IndexedIterator BlasSubmatrix<_Field>::IndexedBegin ()
	{
		return IndexedIterator (_M->Begin () + ( (_off) ),
					   _col , _stride, 0, 0);
	}

	template <class _Field>
	typename BlasSubmatrix<_Field>::IndexedIterator BlasSubmatrix<_Field>::IndexedEnd ()
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
	template <class _Field>
	class BlasSubmatrix<_Field>::ConstIndexedIterator {
	public:
		ConstIndexedIterator (){}

		ConstIndexedIterator (const typename BlasMatrix<_Field>::ConstIterator& cur,
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
			else {
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
			else {
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

		const _Field& operator*() const
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

		const _Field & value() const
		{
			return this->operator*();

		}

	protected:
		typename BlasMatrix<_Field>::ConstIterator _cur;
		size_t _stride;
		size_t _c_dim;
		size_t _r_idx;
		size_t _c_idx;
	};

	template <class _Field>
	typename BlasSubmatrix<_Field>::ConstIndexedIterator BlasSubmatrix<_Field>::IndexedBegin () const
	{
		return ConstIndexedIterator (_M->Begin () + ( _off ),
						_row, _stride, 0, 0);
	}

	template <class _Field>
	typename BlasSubmatrix<_Field>::ConstIndexedIterator BlasSubmatrix<_Field>::IndexedEnd () const
	{
		return ConstIndexedIterator (_M->Begin () + ( (_row) * _stride + (_off+_col) ),
						_col, _stride, _row-1, _col-1);
	}

	////////
	template <class _Field>
	typename BlasSubmatrix<_Field>::RowIterator BlasSubmatrix<_Field>::rowBegin ()
	{
		return RowIterator (_M->Begin () + ( _off  ),
				    _col, _stride);
	}

	template <class _Field>
	typename BlasSubmatrix<_Field>::RowIterator BlasSubmatrix<_Field>::rowEnd ()
	{
		return RowIterator (_M->Begin () + ( (_row) * _stride + _off ),
				    _col, _stride);
	}

	template <class _Field>
	typename BlasSubmatrix<_Field>::ConstRowIterator BlasSubmatrix<_Field>::rowBegin () const
	{
		return ConstRowIterator (_M->Begin () + ( _off ),
					 _col, _stride);
	}

	template <class _Field>
	typename BlasSubmatrix<_Field>::ConstRowIterator BlasSubmatrix<_Field>::rowEnd () const
	{
		return ConstRowIterator (_M->Begin () + ( (_row) * _stride + _off ),
					 _col, _stride);
	}

	template <class _Field>
	typename BlasSubmatrix<_Field>::ColIterator BlasSubmatrix<_Field>::colBegin ()
	{
		return ColIterator (_M->Begin () + ( _off ),
				    _stride, _row);
	}

	template <class _Field>
	typename BlasSubmatrix<_Field>::ColIterator BlasSubmatrix<_Field>::colEnd ()
	{
		return ColIterator (_M->Begin () + ( (_col) + _off ),
				    _stride, _row);
	}

	template <class _Field>
	typename BlasSubmatrix<_Field>::ConstColIterator BlasSubmatrix<_Field>::colBegin () const
	{
		return ConstColIterator (_M->Begin () + ( _off ),
					 _stride, _row);
	}

	template <class _Field>
	typename BlasSubmatrix<_Field>::ConstColIterator BlasSubmatrix<_Field>::colEnd () const
	{
		return ConstColIterator (_M->Begin () + ( (_col) + _off ),
					 _stride, _row);
	}

	/*  operators */
	template <class _Field>
	typename BlasSubmatrix< _Field>::Row BlasSubmatrix< _Field>::operator[] (size_t i)
	{
		return Row (_M.Begin () + (_r0+i) * _stride, _M.Begin () + ((_r0+i) * _stride + _stride) );
	}

	template <class _Field>
	typename BlasSubmatrix< _Field>::ConstRow BlasSubmatrix< _Field>::operator[] (size_t i) const
	{
		return Row (_M.Begin () + (_r0+i) * _stride, _M.Begin () + ((_r0+i) * _stride + _stride) );
	}

} // LinBox

//////////////////
//     MISC     //
//////////////////

namespace LinBox
{
} // LinBox
#endif // __LINBOX_blas_submatrix_INL
