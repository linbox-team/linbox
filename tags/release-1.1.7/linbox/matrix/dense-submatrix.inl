/* linbox/matrix/dense-submatrix.inl
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
 * -----------------------------------------------------------
 * 2002-10-27  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 * 
 * Rename from densesubmatrix.C
 *
 * Constructor modifications: changed the interface to match Submatrix
 * -----------------------------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __LINBOX_dense_submatrix_INL
#define __LINBOX_dense_submatrix_INL

#include "linbox/util/debug.h"
#include "linbox/matrix/dense.h"
#include "linbox/matrix/dense-submatrix.h"

namespace LinBox
{

template <class _Element>
DenseSubmatrix<_Element>::DenseSubmatrix (DenseMatrixBase<_Element> &M,
					 size_t row,
					 size_t col,
					 size_t rowdim,
					 size_t coldim)
	: _M (&M), _beg_row (row), _end_row (row + rowdim), _beg_col (col), _end_col (col + coldim)
{
	linbox_check (_beg_row <= _end_row);
	linbox_check (_beg_col <= _end_col);
	linbox_check (_end_row <= M.rowdim ());
	linbox_check (_end_col <= M.coldim ());
} 


template <class _Element>
DenseSubmatrix<_Element>::DenseSubmatrix (DenseMatrixBase<_Element> &M)
	: _M(&M), _beg_row(0), _end_row(M.rowdim()), _beg_col(0), _end_col(M.coldim()) {}



template <class _Element>
DenseSubmatrix<_Element>::DenseSubmatrix (const DenseSubmatrix<_Element> &SM,
					 size_t row,
					 size_t col,
					 size_t rowdim,
					 size_t coldim)
	: _M (SM._M),
	  _beg_row (SM._beg_row + row),
	  _end_row (SM._beg_row + row + rowdim),
	  _beg_col (SM._beg_col + col),
	  _end_col (SM._beg_col + col + coldim)
{
	linbox_check (_beg_row <= _end_row);
	linbox_check (_beg_col <= _end_col);
	linbox_check (_end_row - _beg_row <= SM.rowdim ());
	linbox_check (_end_col - _beg_col <= SM.coldim ());
}
  
template <class _Element>
DenseSubmatrix<_Element>::DenseSubmatrix (const DenseSubmatrix<_Element> &SM)
	: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col)
{
}

template <class _Element>
DenseSubmatrix<_Element>& DenseSubmatrix<_Element>::operator=(const DenseSubmatrix<_Element> &SM)
{
	_M = (SM._M);
	_beg_row = SM._beg_row;
	_end_row = SM._end_row;
	_beg_col = SM._beg_col;
	_end_col = SM._end_col;

	return *this;
}
  
template <class _Element>
class DenseSubmatrix<_Element>::RawIterator
{
    public:
	RawIterator (){}

	RawIterator (const typename DenseMatrixBase<_Element>::RawIterator& cur,
		     const size_t c_dim,
		     const size_t stride,
		     const size_t c_idx)
		: _cur (cur), _c_dim (c_dim), _stride(stride), _c_idx (c_idx)
	{}

	RawIterator& operator = (const RawIterator& r)
	{
		_cur = r._cur;
		_c_dim = r._c_dim;
		_stride = r._stride;
		_c_idx = r._c_idx;
		return *this;
	}

	RawIterator& operator ++()
	{
		if (_c_idx < _c_dim - 1){
			++_cur; ++_c_idx;
		} else {
			_cur = _cur + _stride - _c_dim + 1;
			_c_idx = 0;
		}

		return *this;  
	}

	RawIterator& operator++ (int)
	{
		return this->operator++ ();
	}

	bool operator != (const RawIterator& r) const
	{
		return (_cur != r._cur || _c_dim != r._c_dim) || (_stride != r._stride) || (_c_idx != r._c_idx);
	} 

	_Element& operator * ()
		{ return *_cur; }

	const _Element& operator * () const
		{ return *_cur; }

    protected:
	typename DenseMatrixBase<_Element>::RawIterator _cur;
	size_t _c_dim;
	size_t _stride;
	size_t _c_idx;
};
  
template <class _Element>
class DenseSubmatrix<_Element>::ConstRawIterator
{
    public:
	ConstRawIterator (){}

	ConstRawIterator (const typename DenseMatrixBase<_Element>::ConstRawIterator& cur,
		     const size_t c_dim,
		     const size_t stride,
		     const size_t c_idx)
		: _cur (cur), _c_dim (c_dim), _stride(stride), _c_idx (c_idx)
	{}

	ConstRawIterator& operator = (const ConstRawIterator& r)
	{
		_cur = r._cur;
		_c_dim = r._c_dim;
		_stride = r._stride;
		_c_idx = r._c_idx;
		return *this;
	}

	ConstRawIterator& operator ++()
	{
		if (_c_idx < _c_dim - 1){
			++_cur; ++_c_idx;
		} else {
			_cur = _cur + _stride - _c_dim + 1;
			_c_idx = 0;
		}

		return *this;  
	}

	ConstRawIterator& operator++ (int)
	{
		return this->operator++ ();
	}

	bool operator != (const ConstRawIterator& r) const
	{
		return (_cur != r._cur) || (_c_dim != r._c_dim) || (_stride != r._stride) || (_c_idx != r._c_idx);
	} 

	const _Element& operator * () const
		{ return *_cur; }

    protected:
	typename DenseMatrixBase<_Element>::ConstRawIterator _cur;
	size_t _c_dim;
	size_t _stride;
	size_t _c_idx;
};
  
// template <class Element>
// class DenseSubmatrix<Element>::ConstRawIterator
// {   
//     public:
// 	ConstRawIterator (){}

// 	ConstRawIterator ( const typename DenseMatrixBase<Element>::ConstRawIterator& cur, 
// 			  size_t cont_len,
// 			  size_t gap_len)
// 		:   _beg (beg), _cur (cur), _cont_len (cont_len), _gap_len (gap_len)
// 	{}

// 	ConstRawIterator& operator = (const RawIterator& r)
// 	{
// 		_cur = r._cur;
// 		_beg = r._beg;
// 		_cont_len = r._cont_len;
// 		_gap_len = r._gap_len;
// 		return *this;
// 	}

// 	ConstRawIterator& operator = (const ConstRawIterator& r)
// 	{
// 		_cur = r._cur;
// 		_beg = r._beg;
// 		_cont_len = r._cont_len;
// 		_gap_len = r._gap_len;
// 		return *this;
// 	}

// 	ConstRawIterator& operator++()
// 	{
// 		if (((_cur - _beg + 1) % _cont_len) != 0)
// 			++_cur;
// 		else
// 		{
// 			_cur = _cur + _gap_len + 1;
// 			_beg = _beg + _gap_len + _cont_len;
// 		}
// 		return *this;
// 	}

// 	ConstRawIterator operator++(int)
// 	{
// 		ConstRawIterator tmp = *this;
// 		this->operator++();
// 		return tmp;
// 	}

// 	bool operator != (const ConstRawIterator& r) const
// 	{
// 		return (_cur != r._cur) || (_beg != r._beg) || (_cont_len != r._cont_len) || (_gap_len != r._gap_len);
// 	}
    
// 	const Element& operator*()
// 	{ return *_cur; }

// // 	Element& operator*()
// // 		{ return *_cur; }

// 	const Element& operator*() const
// 	{ return *_cur; }

//     protected:
// 	typename DenseMatrixBase<Element>::ConstRawIterator _beg;
// 	typename DenseMatrixBase<Element>::ConstRawIterator _cur;
// 	size_t _cont_len;
// 	size_t _gap_len;
// };
  
template <class _Element>
typename DenseSubmatrix<_Element>::RawIterator DenseSubmatrix<_Element>::rawBegin ()
{
	return RawIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col, 
			    coldim (), _M->coldim (), 0);
}
  
template <class _Element>
typename DenseSubmatrix<_Element>::RawIterator DenseSubmatrix<_Element>::rawEnd ()
{
	return RawIterator (_M->rawBegin () + _end_row * _M->coldim () + _beg_col, 
			    coldim (), _M->coldim (), 0);
}
template <class _Element>
typename DenseSubmatrix<_Element>::ConstRawIterator DenseSubmatrix<_Element>::rawBegin () const
{
	return ConstRawIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col, 
			    coldim (), _M->coldim (), 0);
}
  
template <class _Element>
typename DenseSubmatrix<_Element>::ConstRawIterator DenseSubmatrix<_Element>::rawEnd () const 
{
	return ConstRawIterator (_M->rawBegin () + _end_row * _M->coldim () + _beg_col, 
			    coldim (), _M->coldim (), 0);
}
    
// template <class Element>
// typename DenseSubmatrix<Element>::ConstRawIterator DenseSubmatrix<Element>::rawBegin () const
// {
// 	return ConstRawIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col, 
// 				 _M->rawBegin () + _beg_row * _M->coldim () + _beg_col, 
// 				 coldim (), _M->coldim () - coldim ());
// }
  
// template <class Element>
// typename DenseSubmatrix<Element>::ConstRawIterator DenseSubmatrix<Element>::rawEnd () const
// {
// 	return ConstRawIterator (_M->rawBegin () + _end_row * _M->coldim () + _beg_col, 
// 				 _M->rawBegin () + _end_row * _M->coldim () + _beg_col, 
// 				 coldim (), _M->coldim () - coldim ());
// }
//////
template <class _Element>
class DenseSubmatrix<_Element>::RawIndexedIterator
{   
       public:
	RawIndexedIterator (){}
	
	RawIndexedIterator (const typename DenseMatrixBase<_Element>::RawIterator& cur, 
			    size_t c_dim,
			    size_t stride,
			    size_t r_idx,
			    size_t c_idx)
		:   _cur (cur), _c_dim (c_dim), _stride (stride), _r_idx(r_idx), _c_idx (c_idx)
	{}
	
	RawIndexedIterator& operator = (const RawIndexedIterator& r)
	{
		_cur = r._cur;
		_stride = r._stride;
		_c_dim = r._c_dim;
		_r_idx = r._r_idx;
		_c_idx = r._c_idx;
		return *this;
	}
	
	RawIndexedIterator& operator++()
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
	
	RawIndexedIterator& operator--()
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

	RawIndexedIterator operator++(int)
	{
		RawIndexedIterator tmp = *this;
		this->operator++();
		return tmp;
	}

	RawIndexedIterator operator--(int)
	{
		RawIndexedIterator tmp = *this;
		this->operator--();
		return tmp;
	}

	bool operator != (const RawIndexedIterator& r) const
	{
		return ((_c_idx != r._c_idx) || (_r_idx != r._r_idx) ||(_stride != r._stride) || (_c_dim != r._c_dim) );
	}
    
	const _Element& operator*() const {return *_cur;}

	_Element& operator*() {return *_cur;}

	size_t rowIndex () const { return _r_idx; }

	size_t colIndex () const { return _c_idx; }

	const _Element& value () const {return *_cur;}

     protected:
	typename DenseMatrixBase<_Element>::RawIterator _cur;
	size_t _stride;
	size_t _c_dim;
	size_t _r_idx;
	size_t _c_idx;
};
  
template <class _Element>
typename DenseSubmatrix<_Element>::RawIndexedIterator DenseSubmatrix<_Element>::rawIndexedBegin ()
{
	return RawIndexedIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col, 
				   coldim(), _M->coldim(), 0, 0);
}
  
template <class _Element>
typename DenseSubmatrix<_Element>::RawIndexedIterator DenseSubmatrix<_Element>::rawIndexedEnd ()
{
	return RawIndexedIterator (_M->rawBegin () + _end_row * _M->coldim () + _end_col, 
				   coldim(), _M->coldim (), rowdim()-1, coldim()-1);
}

template <class _Element>
class DenseSubmatrix<_Element>::ConstRawIndexedIterator
{   
    public:
	ConstRawIndexedIterator (){}

	ConstRawIndexedIterator (const typename DenseMatrixBase<_Element>::ConstRawIterator& cur, 
				 size_t c_dim,
				 size_t stride,
				 size_t r_idx,
				 size_t c_idx)
		:   _cur (cur), _stride (stride), _c_dim (c_dim), _r_idx(r_idx), _c_idx (c_idx)
	{}

	ConstRawIndexedIterator& operator = (const RawIndexedIterator& r)
	{
		_cur = r._cur;
		_stride = r._stride;
		_c_dim = r._c_dim;
		_r_idx = r._r_idx;
		_c_idx = r._c_idx;
		return *this;
	}

	ConstRawIndexedIterator& operator = (const ConstRawIndexedIterator& r)
	{
		_cur = r._cur;
		_stride = r._stride;
		_c_dim = r._c_dim;
		_r_idx = r._r_idx;
		_c_idx = r._c_idx;
		return *this;
	}

	ConstRawIndexedIterator& operator++()
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
	RawIndexedIterator& operator--()
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

	ConstRawIndexedIterator operator++(int)
	{
		ConstRawIndexedIterator tmp = *this;
		this->operator++();
		return tmp;
	}

	ConstRawIndexedIterator operator--(int)
	{
		ConstRawIndexedIterator tmp = *this;
		this->operator--();
		return tmp;
	}

	size_t rowIndex () const { return _r_idx; }

	size_t colIndex () const { return _c_idx; }

	bool operator != (const ConstRawIndexedIterator& r) const
	{
		return ((_c_idx != r._c_idx) || (_r_idx != r._r_idx) ||(_stride != r._stride) || (_c_dim != r._c_dim) );
	}
	
	const _Element& operator*() const {return *_cur;}


    	friend std::ostream& operator<<(std::ostream& out, const ConstRawIndexedIterator m) {
            return out /* << m._cur << ' ' */
                       << m._stride << ' '
                       << m._c_dim << ' '
                       << m._r_idx << ' ' 
                       << m._c_idx;
        }    

    protected:
	typename DenseMatrixBase<_Element>::ConstRawIterator _cur;
	size_t _stride;
	size_t _c_dim;
	size_t _r_idx;
	size_t _c_idx;
};

template <class _Element>
typename DenseSubmatrix<_Element>::ConstRawIndexedIterator DenseSubmatrix<_Element>::rawIndexedBegin () const
{
	return ConstRawIndexedIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col, 
					coldim(), _M->coldim(), 0, 0);
}
  
template <class _Element>
typename DenseSubmatrix<_Element>::ConstRawIndexedIterator DenseSubmatrix<_Element>::rawIndexedEnd () const
{
	return ConstRawIndexedIterator (_M->rawBegin () + _end_row * _M->coldim () + _end_col, 
					coldim (), _M->coldim (), rowdim()-1, coldim()-1);
}

////////  
template <class _Element>
typename DenseSubmatrix<_Element>::RowIterator DenseSubmatrix<_Element>::rowBegin ()
{
	return RowIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col,
			    _end_col - _beg_col, _M->coldim ());
}
 
template <class _Element>
typename DenseSubmatrix<_Element>::RowIterator DenseSubmatrix<_Element>::rowEnd ()
{
	return RowIterator (_M->rawBegin () + _end_row * _M->coldim () + _beg_col,
			    _end_col - _beg_col, _M->coldim ());
}

template <class _Element>
typename DenseSubmatrix<_Element>::ConstRowIterator DenseSubmatrix<_Element>::rowBegin () const
{
	return ConstRowIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col,
				 _end_col - _beg_col, _M->coldim ());
}
  
template <class _Element>
typename DenseSubmatrix<_Element>::ConstRowIterator DenseSubmatrix<_Element>::rowEnd () const
{
	return ConstRowIterator (_M->rawBegin () + _end_row * _M->coldim () + _beg_col,
				 _end_col - _beg_col, _M->coldim ());
}

template <class _Element>
typename DenseSubmatrix<_Element>::ColIterator DenseSubmatrix<_Element>::colBegin ()
{
	return ColIterator (_M->rawBegin () + _beg_col + _beg_row * _M->coldim (),
			    _M->coldim (), rowdim ());
}

template <class _Element>
typename DenseSubmatrix<_Element>::ColIterator DenseSubmatrix<_Element>::colEnd ()
{
	return ColIterator (_M->rawBegin () + _end_col + _beg_row * _M->coldim (),
			    _M->coldim (), rowdim ());
}

template <class _Element>
typename DenseSubmatrix<_Element>::ConstColIterator DenseSubmatrix<_Element>::colBegin () const
{
	return ConstColIterator (_M->rawBegin () + _beg_col + _beg_row * _M->coldim (),
				 _M->coldim (), rowdim ());
}

template <class _Element>
typename DenseSubmatrix<_Element>::ConstColIterator DenseSubmatrix<_Element>::colEnd () const
{
	return ConstColIterator (_M->rawBegin () + _end_col + _beg_row * _M->coldim (),
				 _M->coldim (), rowdim ());
}

template <class _Element>
template <class Field>
std::istream& DenseSubmatrix<_Element>::read (std::istream &file, const Field& field)
{
	RawIterator p;

	for (p = rawBegin (); p != rawEnd (); ++p) {
		// each entry is seperated by one space.
		file.ignore (1);
		field.read (file, *p);
	}

	return file;
}

template <class _Element>
template <class Field>
std::ostream &DenseSubmatrix<_Element>::write (std::ostream &os, const Field& field, bool mapleFormat) const
{
	ConstRowIterator p;

	integer c;
	//int wid;

	field.cardinality (c);
	//wid = (int) ceil (log ((double) c) / M_LN10); //BB : not used !

	typename ConstRow::const_iterator pe;

	if (mapleFormat) os << "[";

	for (p = rowBegin (); p != rowEnd (); ++p) {
		if (mapleFormat && (p != rowBegin()))
			os << ',';
		if (mapleFormat) os << "[";

		for (pe = p->begin (); pe != p->end (); ++pe) {
			if (mapleFormat && (pe != p->begin())) os << ',';
			// matrix base does not provide this field(), maybe should?
			//_M.field ().write (os, *pe);
		        //os << *pe;
			//fixed by using extra field

			field.write (os, *pe);
			os << " ";
		}

		if (!mapleFormat)
			os << std::endl;
		else os << ']';
	}

	if (mapleFormat) os << ']';
	return os;
}

} // namespace LinBox

#endif // __LINBOX_dense_submatrix_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
