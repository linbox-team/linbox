/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/dense-base.inl
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
 *
 * See COPYING for license information
 */

#ifndef __DENSE_BASE_INL
#define __DENSE_BASE_INL

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#include "linbox/blackbox/dense-base.h"
#include "linbox/util/debug.h"

namespace LinBox
{

template <class Element>
class DenseMatrixBase<Element>::ConstColOfRowsIterator
{
    public:
	ConstColOfRowsIterator (const typename Rep::const_iterator& p, size_t len, size_t d)
		: _row (p, p + len), _dis (d) {}
    
	ConstColOfRowsIterator () {}
    
	ConstColOfRowsIterator (const ConstColOfRowsIterator& colp)
		: _row (colp._row), _dis (colp._dis) {}
    
	ConstColOfRowsIterator& operator = (const ConstColOfRowsIterator& colp)
	{
		_row = colp._row;
		_dis = colp._dis;
		return *this;
	}

	ConstColOfRowsIterator& operator++ ()
	{
		_row = ConstRow (_row.begin () + _dis, _row.end () + _dis);
		return *this;
	}

	ConstColOfRowsIterator  operator++ (int)
	{
		ColOfRowsIterator tmp (*this);
		++_row;
		return tmp;
	}

	ConstColOfRowsIterator& operator+ (int i)
	{
		_row = ConstRow (_row.begin () + _dis * i, _row.end () + _dis * i);
		return *this;
	}

	ConstRow operator[] (int i) const
		{ return ConstRow (_row.begin () + _dis * i, _row.end () + _dis * i); }

	ConstRow* operator-> ()
		{ return &_row; }

	ConstRow& operator* ()
		{ return _row; }
    
	bool operator!= (const ConstColOfRowsIterator& c) const
		{ return (_row.begin () != c._row.begin ()) || (_row.end () != c._row.end ()) || (_dis != c._dis); }

    private:
	ConstRow _row;
	size_t _dis;
};

template <class Element>
class DenseMatrixBase<Element>::ColOfRowsIterator
{
    public:
	ColOfRowsIterator (const typename Rep::iterator& p, size_t len, size_t d)
		: _row (p, p + len), _dis (d){}

	ColOfRowsIterator () {}

	ColOfRowsIterator (const ColOfRowsIterator& colp)
		: _row (colp._row), _dis (colp._dis) {}
    
	ColOfRowsIterator& operator= (const ColOfRowsIterator& colp)
	{
		_row = colp._row;
		_dis = colp._dis;
		return *this;
	}
    
    
	ColOfRowsIterator& operator++ ()
	{
		_row = Row (_row.begin () + _dis, _row.end () + _dis);
		return *this;
	}
    
	ColOfRowsIterator  operator++ (int)
	{
		ColOfRowsIterator tmp (*this);
		++_row;
		return tmp;
	}
    
	ColOfRowsIterator& operator+ (int i)
	{
		_row = tRow (_row.begin () + _dis * i, _row.end () + _dis * i);
		return *this;
	}

	Row operator[] (int i) const
		{ return Row (const_cast<Row&> (_row).begin () + _dis * i,
			      const_cast<Row&> (_row).end () + _dis * i); }

	Row* operator-> ()
		{ return &_row; }
    
	Row& operator* ()
		{ return _row; }
 
	bool operator!= (const ColOfRowsIterator& c) const
		{ return (_row.begin () != c._row.begin ()) || (_row.end () != c._row.end ()) || (_dis != c._dis); }

	operator ConstColOfRowsIterator ()
		{ return ConstColOfRowsIterator (_row.begin (), _row.size (), _dis); }

    private:
	Row _row;
	size_t _dis;
};

template <class Element>
class DenseMatrixBase<Element>::ConstRowOfColsIterator
{
    public:
	ConstRowOfColsIterator (typename Rep::const_iterator p, size_t stride, size_t len)
		: _col (Subiterator<typename Rep::const_iterator> (p, stride),
			Subiterator<typename Rep::const_iterator> (p + len * stride, stride)), _stride (stride)
	{}
    
	ConstRowOfColsIterator (const ConstCol& col, size_t stride)
		:_col (col), _stride (stride){}

	ConstRowOfColsIterator () {}
    
	ConstRowOfColsIterator (const ConstRowOfColsIterator& rowp)
		:_col (rowp._col){}

	ConstRowOfColsIterator& operator= (const ConstRowOfColsIterator& rowp)
	{
		_col = rowp._col;
		_stride = rowp._stride;
		return *this;
	}

	ConstRowOfColsIterator& operator++ ()
	{
		_col = ConstCol (Subiterator<typename Rep::const_iterator> (_col.begin ().operator-> () + 1, _stride),
				 Subiterator<typename Rep::const_iterator> (_col.end ().operator-> () + 1, _stride));
		return *this;
	}

	ConstRowOfColsIterator  operator++ (int)
	{
		Col tmp (_col);
		this->operator++ ();
		return tmp;
	}

	ConstRowOfColsIterator& operator+ (int i)
	{ 
		_col = ConstCol (Subiterator<typename Rep::const_iterator> (_col.begin ().operator-> () + i, _stride),
				 Subiterator<typename Rep::const_iterator> (_col.end ().operator-> () + i, _stride));
		return *this;
	}


	ConstCol operator[] (int i) const
		{ return ConstCol (Subiterator<typename Rep::const_iterator> (_col.begin ().operator-> () + i, _stride), 
				   Subiterator<typename Rep::const_iterator> (_col.end ().operator-> () + i, _stride)); }

	ConstCol* operator-> ()
		{ return &_col; }
 
	ConstCol& operator* ()
		{ return _col; }
    
	bool operator!= (const ConstRowOfColsIterator& c) const
		{ return (_col.begin () != c._col.begin ()) || (_col.end () != c._col.end ()); }
    
    private:
	ConstCol _col;
	size_t _stride;
};

template <class Element>
class DenseMatrixBase<Element>::RowOfColsIterator
{
    public:
	RowOfColsIterator (typename Rep::iterator p, size_t stride, size_t len)
		: _col (Subiterator<typename Rep::iterator> (p, stride),
			Subiterator<typename Rep::iterator> (p+len*stride, stride)), _stride (stride)
	{}
    
	RowOfColsIterator () {}
    
	RowOfColsIterator (const RowOfColsIterator& rowp)
		:_col (rowp._col){}
    
	RowOfColsIterator& operator= (const RowOfColsIterator& rowp)
	{
		_col = rowp._col;
		_stride = rowp._stride;
		return *this;
	}
    
	const RowOfColsIterator& operator= (const RowOfColsIterator& rowp) const
	{
		const_cast<RowOfColsIterator*> (this)->_col = rowp._col;
		return *this;
	}
    
	RowOfColsIterator& operator++ ()
	{
		_col = Col (Subiterator<typename Rep::iterator> (_col.begin ().operator-> () + 1, _stride),
			    Subiterator<typename Rep::iterator> (_col.end ().operator-> () + 1, _stride));
		return *this;
	}
    
	RowOfColsIterator  operator++ (int)
	{
		Col tmp (_col);
		this->operator++ ();
		return tmp;
	}
        
	RowOfColsIterator& operator+ (int i)
	{ 
		_col = Col (Subiterator<typename Rep::const_iterator> (_col.begin ().operator-> () + i, _stride), 
			    Subiterator<typename Rep::const_iterator> (_col.end ().operator-> () + i, _stride));
		return *this;
	}
    

	Col operator[] (int i) const
		{ return Col (Subiterator<typename Rep::iterator> (const_cast<Col&> (_col).begin ().operator-> () + i, _stride), 
			      Subiterator<typename Rep::iterator> (const_cast<Col&> (_col).end ().operator-> () + i, _stride)); }
    
	Col* operator-> ()
		{ return &_col; }
    
	Col& operator* ()
		{ return _col; }
    
	bool operator!= (const RowOfColsIterator& c) const
		{ return (_col.begin () != c._col.begin ()) || (_col.end () != c._col.end ()); }
    
	operator ConstRowOfColsIterator ()
	{
		ConstCol tmp;
		tmp = _col;
		return ConstRowOfColsIterator (tmp, _stride);
	}
    
    private:
	Col _col;
	size_t _stride;
};

/// entry access raw view.  Size m*n vector in C (row major) order.
template <class Element>
typename DenseMatrixBase<Element>::RawIterator DenseMatrixBase<Element>::rawBegin ()
	{ return _rep.begin (); }  

template <class Element>
typename DenseMatrixBase<Element>::RawIterator DenseMatrixBase<Element>::rawEnd ()
	{ return _rep.end (); }
  
template <class Element>
typename DenseMatrixBase<Element>::ConstRawIterator DenseMatrixBase<Element>::rawBegin () const
	{ return _rep.begin (); }  

template <class Element>
typename DenseMatrixBase<Element>::ConstRawIterator DenseMatrixBase<Element>::rawEnd () const
	{ return _rep.end (); }

template <class Element>
typename DenseMatrixBase<Element>::ColOfRowsIterator DenseMatrixBase<Element>::colOfRowsBegin ()
	{ return ColOfRowsIterator (_rep.begin (), _cols, _cols); }

template <class Element>
typename DenseMatrixBase<Element>::ColOfRowsIterator DenseMatrixBase<Element>::colOfRowsEnd ()
	{ return ColOfRowsIterator (_rep.end (), _cols, _cols); }
  
template <class Element>
typename DenseMatrixBase<Element>::ConstColOfRowsIterator DenseMatrixBase<Element>::colOfRowsBegin () const
	{ return ConstColOfRowsIterator (_rep.begin (), _cols, _cols); }  

template <class Element>
typename DenseMatrixBase<Element>::ConstColOfRowsIterator DenseMatrixBase<Element>::colOfRowsEnd () const
	{return ConstColOfRowsIterator (_rep.end (), _cols, _cols); }
  
template <class Element>
typename DenseMatrixBase<Element>::RowOfColsIterator DenseMatrixBase<Element>::rowOfColsBegin ()
	{ return  DenseMatrixBase<Element>::RowOfColsIterator (_rep.begin (), _cols, _rows); }

template <class Element>
typename DenseMatrixBase<Element>::RowOfColsIterator DenseMatrixBase<Element>::rowOfColsEnd ()
	{ return  DenseMatrixBase<Element>::RowOfColsIterator (_rep.begin ()+_cols, _cols, _rows); }
  
template <class Element>
typename DenseMatrixBase<Element>::ConstRowOfColsIterator DenseMatrixBase<Element>::rowOfColsBegin () const
	{ return  DenseMatrixBase<Element>::ConstRowOfColsIterator (_rep.begin (), _cols, _rows); }
  
template <class Element>
typename DenseMatrixBase<Element>::ConstRowOfColsIterator DenseMatrixBase<Element>::rowOfColsEnd () const
	{ return  DenseMatrixBase<Element>::ConstRowOfColsIterator (_rep.begin ()+_cols, _cols, _rows); }
  
template <class Element>
template <class Field>
void DenseMatrixBase<Element>::read (const Field &F, std::istream &file)
{
	RawIterator p;

	for (p = rawBegin (); p != rawEnd (); ++p) {
		file.ignore (1);
		F.read (file, *p);
	}
}
  
template <class Element>
template <class Field>
std::ostream& DenseMatrixBase<Element>::write (const Field &F, std::ostream &os) const
{
	ConstColOfRowsIterator p;

	integer c;
	int wid;

	F.cardinality (c);
	wid = (int) ceil (log ((double) c) / M_LN10);

	for (p = colOfRowsBegin (); p != colOfRowsEnd (); ++p) {
		typename ConstRow::const_iterator pe;

		commentator.indent (os);

		os << "  [ ";

		for (pe = p->begin (); pe != p->end (); ++pe) {
			os.width (wid);
			F.write (os, *pe);
			os << ' ';
		}

		os << " ]" << endl;
	}

	return os;
}

template <class Element>
class DenseMatrixBase<Element>::RawIndexedIterator
{
	mutable size_t _r_index;
	mutable size_t _c_index;
	mutable size_t _dim;
	typename Rep::iterator _begin;

    public:
	RawIndexedIterator (const size_t  &dim,
			  const size_t  &r_index,
			  const size_t  &c_index,
			  const typename Rep::iterator &begin)
		: _r_index (r_index), _c_index (c_index), _dim (dim), _begin (begin)
	{}
	
	RawIndexedIterator ():_r_index (0), _c_index (0), _dim (1), _begin (0){}

	RawIndexedIterator (const RawIndexedIterator& r)
		: _r_index (r._r_index), _c_index (r._c_index), _dim (r._dim), _begin (r._begin)
	{}

	RawIndexedIterator& operator = (const RawIndexedIterator &iter)
	{
		_r_index = iter._r_index;
		_c_index = iter._c_index;
		_dim = iter._dim;
		_begin = iter._begin;
		return *this;
	}
	
	bool operator == (const RawIndexedIterator &iter) const
		{ return (_r_index == iter._r_index) &&
			  (_c_index == iter._c_index) &&
			  (_dim == iter._dim) &&
			  (_begin==iter._begin); }

	bool operator != (const RawIndexedIterator& iter) const
		{ return (_r_index != iter._r_index) ||
			  (_c_index != iter._c_index) ||
			  (_dim != iter._dim) ||
			  (_begin!=iter._begin); }
	
	RawIndexedIterator &operator ++ ()
	{
		++_c_index;

		if (_c_index == _dim) {
			_c_index = 0;
			++_r_index;
		}

		return *this;
	}

	const RawIndexedIterator &operator ++ () const
	{
		++_c_index;

		if (_c_index==_dim) {
			_c_index = 0;
			++_r_index;
		}

		return *this;
	}
	
	RawIndexedIterator operator ++ (int) const
	{
		RawIndexedIterator tmp = *this;
		++(*this);
		return tmp;
	}

	RawIndexedIterator &operator -- ()
	{ 
		if (_c_index)
			--_c_index;
		else {
			--_r_index;
			_c_index = _dim - 1;
		}

		return *this;
	}

	const RawIndexedIterator &operator -- () const
	{ 
		if (_c_index)
			--_c_index;
		else {
			--_r_index;
			_c_index = _dim - 1;
		}

		return *this;
	}

	RawIndexedIterator operator -- (int) const
	{
		RawIndexedIterator tmp = *this;
		--(*this);
		return tmp;
	}	

	Element &operator * ()
		{ return *(_begin + (_r_index * _dim + _c_index)); }

	const Element &operator * () const
		{ return *(_begin + (_r_index * _dim + _c_index)); }
	
	Element* operator -> ()
		{ return _begin + (_r_index * _dim + _c_index); }
	
	const Element *operator -> () const
		{ return _begin + (_r_index * _dim + _c_index); }

	size_t rowIndex () const
		{ return _r_index; }

	size_t colIndex () const
		{ return _c_index; }
};

template <class Element>
typename DenseMatrixBase<Element>::RawIndexedIterator DenseMatrixBase<Element>::rawIndexedBegin () 
{
	return RawIndexedIterator (coldim (), 0, 0, _rep.begin ());
}

template <class Element>
typename DenseMatrixBase<Element>::RawIndexedIterator DenseMatrixBase<Element>::rawIndexedEnd ()
{
	return RawIndexedIterator (coldim (), rowdim (), 0, _rep.begin ());
}
template <class Element>
typename DenseMatrixBase<Element>::ConstRawIndexedIterator DenseMatrixBase<Element>::rawIndexedBegin () const
{
	return RawIndexedIterator (coldim (), 0, 0, _rep.begin ());
}

template <class Element>
typename DenseMatrixBase<Element>::ConstRawIndexedIterator DenseMatrixBase<Element>::rawIndexedEnd () const
{
	return RawIndexedIterator (coldim (), rowdim (), 0, _rep.begin ());
}

} // namespace LinBox

#endif // __DENSE_BASE_INL
