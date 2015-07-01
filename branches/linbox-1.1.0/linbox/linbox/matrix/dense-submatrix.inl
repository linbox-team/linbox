/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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

#ifndef __DENSE_SUBMATRIX_INL
#define __DENSE_SUBMATRIX_INL

#include "linbox/util/debug.h"
#include "linbox/matrix/dense.h"
#include "linbox/matrix/dense-submatrix.h"

namespace LinBox
{

template <class Element>
DenseSubmatrix<Element>::DenseSubmatrix (DenseMatrixBase<Element> &M,
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


template <class Element>
DenseSubmatrix<Element>::DenseSubmatrix (DenseMatrixBase<Element> &M)
	: _M(&M), _beg_row(0), _end_row(M.rowdim()), _beg_col(0), _end_col(M.coldim()) {}



template <class Element>
DenseSubmatrix<Element>::DenseSubmatrix (const DenseSubmatrix<Element> &SM,
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
  
template <class Element>
DenseSubmatrix<Element>::DenseSubmatrix (const DenseSubmatrix<Element> &SM)
	: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col)
{}

template <class Element>
DenseSubmatrix<Element>& DenseSubmatrix<Element>::operator=(const DenseSubmatrix<Element> &SM)
{
	_M = (SM._M);
	_beg_row = SM._beg_row;
	_end_row = SM._end_row;
	_beg_col = SM._beg_col;
	_end_col = SM._end_col;

	return *this;
}
  
template <class Element>
class DenseSubmatrix<Element>::RawIterator
{
    public:
	RawIterator (){}

	RawIterator (const typename DenseMatrixBase<Element>::RawIterator& beg,
		     const typename DenseMatrixBase<Element>::RawIterator& cur,
		     size_t cont_len,
		     size_t gap_len)
		: _beg (beg), _cur (cur), _cont_len (cont_len), _gap_len (gap_len)
	{}

	RawIterator& operator = (const RawIterator& r)
	{
		_cur = r._cur;
		_beg = r._beg;
		_cont_len = r._cont_len;
		_gap_len = r._gap_len;
		return *this;
	}

	RawIterator& operator ++()
	{
		if (((_cur - _beg + 1) % _cont_len) != 0)
			++_cur;
		else {
			_cur = _cur + _gap_len + 1;
			_beg = _beg + _gap_len + _cont_len;
		}

		return *this;  
	}

	RawIterator& operator++ (int)
	{
		return this->operator++ ();
	}

	bool operator != (const RawIterator& r) const
	{
		return (_cur != r._cur) || (_beg != r._beg) || (_cont_len != r._cont_len) || (_gap_len != r._gap_len);
	} 

	Element& operator * ()
		{ return *_cur; }

	const Element& operator * () const
		{ return *_cur; }

    protected:
	typename DenseMatrixBase<Element>::RawIterator _beg;
	typename DenseMatrixBase<Element>::RawIterator _cur;
	size_t _cont_len;
	size_t _gap_len;
};
  
template <class Element>
class DenseSubmatrix<Element>::ConstRawIterator
{   
    public:
	ConstRawIterator (){}

	ConstRawIterator (const typename DenseMatrixBase<Element>::ConstRawIterator& beg, 
			  const typename DenseMatrixBase<Element>::ConstRawIterator& cur, 
			  size_t cont_len,
			  size_t gap_len)
		:   _beg (beg), _cur (cur), _cont_len (cont_len), _gap_len (gap_len)
	{}

	ConstRawIterator& operator = (const RawIterator& r)
	{
		_cur = r._cur;
		_beg = r._beg;
		_cont_len = r._cont_len;
		_gap_len = r._gap_len;
		return *this;
	}

	ConstRawIterator& operator = (const ConstRawIterator& r)
	{
		_cur = r._cur;
		_beg = r._beg;
		_cont_len = r._cont_len;
		_gap_len = r._gap_len;
		return *this;
	}

	ConstRawIterator& operator++()
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

	ConstRawIterator operator++(int)
	{
		ConstRawIterator tmp = *this;
		this->operator++();
		return tmp;
	}

	bool operator != (const ConstRawIterator& r) const
	{
		return (_cur != r._cur) || (_beg != r._beg) || (_cont_len != r._cont_len) || (_gap_len != r._gap_len);
	}
    
	const Element& operator*()
	{ return *_cur; }

// 	Element& operator*()
// 		{ return *_cur; }

	const Element& operator*() const
	{ return *_cur; }

    protected:
	typename DenseMatrixBase<Element>::ConstRawIterator _beg;
	typename DenseMatrixBase<Element>::ConstRawIterator _cur;
	size_t _cont_len;
	size_t _gap_len;
};
  
template <class Element>
typename DenseSubmatrix<Element>::RawIterator DenseSubmatrix<Element>::rawBegin ()
{
	return RawIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col, 
			    _M->rawBegin () + _beg_row * _M->coldim () + _beg_col, 
			    coldim (), _M->coldim () - coldim ());
}
  
template <class Element>
typename DenseSubmatrix<Element>::RawIterator DenseSubmatrix<Element>::rawEnd ()
{
	return RawIterator (_M->rawBegin () + _end_row * _M->coldim () + _beg_col, 
			    _M->rawBegin () + _end_row * _M->coldim () + _beg_col, 
			    coldim (), _M->coldim () - coldim ());
}
    
template <class Element>
typename DenseSubmatrix<Element>::ConstRawIterator DenseSubmatrix<Element>::rawBegin () const
{
	return ConstRawIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col, 
				 _M->rawBegin () + _beg_row * _M->coldim () + _beg_col, 
				 coldim (), _M->coldim () - coldim ());
}
  
template <class Element>
typename DenseSubmatrix<Element>::ConstRawIterator DenseSubmatrix<Element>::rawEnd () const
{
	return ConstRawIterator (_M->rawBegin () + _end_row * _M->coldim () + _beg_col, 
				 _M->rawBegin () + _end_row * _M->coldim () + _beg_col, 
				 coldim (), _M->coldim () - coldim ());
}
  
template <class Element>
typename DenseSubmatrix<Element>::RowIterator DenseSubmatrix<Element>::rowBegin ()
{
	return RowIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col,
			    _end_col - _beg_col, _M->coldim ());
}
 
template <class Element>
typename DenseSubmatrix<Element>::RowIterator DenseSubmatrix<Element>::rowEnd ()
{
	return RowIterator (_M->rawBegin () + _end_row * _M->coldim () + _beg_col,
			    _end_col - _beg_col, _M->coldim ());
}

template <class Element>
typename DenseSubmatrix<Element>::ConstRowIterator DenseSubmatrix<Element>::rowBegin () const
{
	return ConstRowIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col,
				 _end_col - _beg_col, _M->coldim ());
}
  
template <class Element>
typename DenseSubmatrix<Element>::ConstRowIterator DenseSubmatrix<Element>::rowEnd () const
{
	return ConstRowIterator (_M->rawBegin () + _end_row * _M->coldim () + _beg_col,
				 _end_col - _beg_col, _M->coldim ());
}

template <class Element>
typename DenseSubmatrix<Element>::ColIterator DenseSubmatrix<Element>::colBegin ()
{
	return ColIterator (_M->rawBegin () + _beg_col + _beg_row * _M->coldim (),
			    _M->coldim (), rowdim ());
}

template <class Element>
typename DenseSubmatrix<Element>::ColIterator DenseSubmatrix<Element>::colEnd ()
{
	return ColIterator (_M->rawBegin () + _end_col + _beg_row * _M->coldim (),
			    _M->coldim (), rowdim ());
}

template <class Element>
typename DenseSubmatrix<Element>::ConstColIterator DenseSubmatrix<Element>::colBegin () const
{
	return ConstColIterator (_M->rawBegin () + _beg_col + _beg_row * _M->coldim (),
				 _M->coldim (), rowdim ());
}

template <class Element>
typename DenseSubmatrix<Element>::ConstColIterator DenseSubmatrix<Element>::colEnd () const
{
	return ConstColIterator (_M->rawBegin () + _end_col + _beg_row * _M->coldim (),
				 _M->coldim (), rowdim ());
}

template <class Element>
template <class Field>
std::istream& DenseSubmatrix<Element>::read (std::istream &file, const Field& field)
{
	RawIterator p;

	for (p = rawBegin (); p != rawEnd (); ++p) {
		// each entry is seperated by one space.
		file.ignore (1);
		field.read (file, *p);
	}

	return file;
}

template <class Element>
template <class Field>
std::ostream &DenseSubmatrix<Element>::write (std::ostream &os, const Field& field, bool mapleFormat) const
{
	ConstRowIterator p;

	integer c;
	int wid;

	field.cardinality (c);
	wid = (int) ceil (log ((double) c) / M_LN10);

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
	os << std::endl;

	return os;
}

} // namespace LinBox

#endif // __DENSE_SUBMATRIX_INL
