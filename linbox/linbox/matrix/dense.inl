/* linbox/matrix/dense.inl
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
 * 2003-01-11  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Move from blackbox/dense-base.inl to matrix/dense.inl
 * --------------------------------------------------------
 * 2002-10-27  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Split out container/iterator functionality into DenseMatrixBase
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __LINBOX_matrix_dense_INL
#define __LINBOX_matrix_dense_INL

#include <iostream>
#include <cmath>

#include "linbox/matrix/dense.h"
#include "linbox/util/debug.h"

namespace LinBox
{

template <class Element>
class DenseMatrixBase<Element>::ConstRowIterator
{
    public:
	ConstRowIterator (const typename Rep::const_iterator& p, size_t len, size_t d)
		: _row (p, p + len), _dis (d) {}
    
	ConstRowIterator () {}
    
	ConstRowIterator (const ConstRowIterator& colp)
		: _row (colp._row), _dis (colp._dis) {}
    
	ConstRowIterator& operator = (const ConstRowIterator& colp)
	{
		_row = colp._row;
		_dis = colp._dis;
		return *this;
	}

	ConstRowIterator& operator --()
	{
		_row = ConstRow (_row.begin () - _dis, _row.end () - _dis);
		return *this;
	}

	ConstRowIterator  operator-- (int)
        {
                ConstRowIterator tmp (*this);
                --*this;
                return tmp;
	}

	
	ConstRowIterator& operator++ ()
	{
		_row = ConstRow (_row.begin () + _dis, _row.end () + _dis);
		return *this;
	}

	ConstRowIterator  operator++ (int)
	{
		ConstRowIterator tmp (*this);
		++*this;
		return tmp;
	}

	ConstRowIterator operator+ (int i)
		{ return ConstRowIterator (_row.begin () + _dis * i, _row.size (), _dis); }

	ConstRowIterator& operator += (int i)
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
    
	bool operator!= (const ConstRowIterator& c) const
		{ return (_row.begin () != c._row.begin ()) || (_row.end () != c._row.end ()) || (_dis != c._dis); }

    private:
	ConstRow _row;
	size_t _dis;
};

template <class Element>
class DenseMatrixBase<Element>::RowIterator
{
    public:
	RowIterator (const typename Rep::iterator& p, size_t len, size_t d)
		: _row (p, p + len), _dis (d){}

	RowIterator () {}

	RowIterator (const RowIterator& colp)
		: _row (colp._row), _dis (colp._dis) {}
    
	RowIterator& operator = (const RowIterator& colp)
	{
		_row = colp._row;
		_dis = colp._dis;
		return *this;
	}
    
	RowIterator& operator ++ ()
	{
		_row = Row (_row.begin () + _dis, _row.end () + _dis);
		return *this;
	}
    
	RowIterator  operator ++ (int)
	{
		RowIterator tmp (*this);
		++*this;
		return tmp;
	}
    
        RowIterator& operator -- ()
        {
                _row = Row (_row.begin () - _dis, _row.end () - _dis);
                return *this;
        }

        RowIterator  operator -- (int)
        {
                RowIterator tmp (*this);
                --*this;
                return tmp;
        }

	RowIterator operator + (int i)
		{ return RowIterator (_row.begin () + _dis * i, _row.size (), _dis); }

	RowIterator& operator += (int i)
	{
		_row = Row (_row.begin () + _dis * i, _row.end () + _dis * i);
		return *this;
	}

	Row operator[] (int i) const
		{ return Row (const_cast<Row&> (_row).begin () + _dis * i,
			      const_cast<Row&> (_row).end () + _dis * i); }

	Row* operator-> ()
		{ return &_row; }
    
	Row& operator* ()
		{ return _row; }
 
	bool operator!= (const RowIterator& c) const
		{ return (_row.begin () != c._row.begin ()) || (_row.end () != c._row.end ()) || (_dis != c._dis); }

	operator ConstRowIterator ()
		{ return ConstRowIterator (_row.begin (), _row.size (), _dis); }

    private:
	Row _row;
	size_t _dis;
};

template <class Element>
class DenseMatrixBase<Element>::ConstColIterator
{
    public:
	ConstColIterator (typename Rep::const_iterator p, size_t stride, size_t len)
		: _col (Subiterator<typename Rep::const_iterator> (p, stride),
			Subiterator<typename Rep::const_iterator> (p + len * stride, stride)), _stride (stride)
	{}
    
	ConstColIterator (const ConstCol& col, size_t stride)
		:_col (col), _stride (stride){}

	ConstColIterator () {}
    
	ConstColIterator (const ConstColIterator& rowp)
		:_col (rowp._col){}

	ConstColIterator& operator= (const ConstColIterator& rowp)
	{
		_col = rowp._col;
		_stride = rowp._stride;
		return *this;
	}

	ConstColIterator& operator++ ()
	{
		_col = ConstCol (Subiterator<typename Rep::const_iterator> (_col.begin ().operator-> () + 1, _stride),
				 Subiterator<typename Rep::const_iterator> (_col.end ().operator-> () + 1, _stride));
		return *this;
	}

	ConstColIterator  operator++ (int)
	{
		ConstColIterator old(*this);
		this->operator++ ();
		return old;
	}

	ConstColIterator operator + (int i)
		{ return ConstColIterator (_col.begin ().operator-> () + i, _stride, _col.size ()); }

	ConstColIterator& operator += (int i)
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
    
	bool operator!= (const ConstColIterator& c) const
		{ return (_col.begin () != c._col.begin ()) || (_col.end () != c._col.end ()); }
    
    private:
	ConstCol _col;
	size_t _stride;
};

template <class Element>
class DenseMatrixBase<Element>::ColIterator
{
    public:
	ColIterator (typename Rep::iterator p, size_t stride, size_t len)
		: _col (Subiterator<typename Rep::iterator> (p, stride),
			Subiterator<typename Rep::iterator> (p + len * stride, stride)), _stride (stride)
	{}
    
	ColIterator () {}
    
	ColIterator (const ColIterator& rowp)
		:_col (rowp._col){}
    
	ColIterator& operator= (const ColIterator& rowp)
	{
		_col = rowp._col;
		_stride = rowp._stride;
		return *this;
	}
    
	const ColIterator& operator= (const ColIterator& rowp) const
	{
		const_cast<ColIterator*> (this)->_col = rowp._col;
		return *this;
	}
    
	ColIterator& operator++ ()
	{
		_col = Col (Subiterator<typename Rep::iterator> (_col.begin ().operator-> () + 1, _stride),
			    Subiterator<typename Rep::iterator> (_col.end ().operator-> () + 1, _stride));
		return *this;
	}
    
	ColIterator  operator++ (int)
	{
		Col tmp (_col);
		this->operator++ ();
		return tmp;
	}
        
	ColIterator operator + (int i)
		{ return ColIterator (_col.begin ().operator-> () + i, _stride, _col.size ()); }

	ColIterator& operator += (int i)
	{
		_col = Col (Subiterator<typename Rep::iterator> (_col.begin ().operator-> () + i, _stride),
			    Subiterator<typename Rep::iterator> (_col.end ().operator-> () + i, _stride));
		return *this;
	}

	Col operator[] (int i) const
		{ return Col (Subiterator<typename Rep::iterator> (const_cast<Col&> (_col).begin ().operator-> () + i, _stride), 
			      Subiterator<typename Rep::iterator> (const_cast<Col&> (_col).end ().operator-> () + i, _stride)); }

	Col* operator-> ()
		{ return &_col; }

	Col& operator* ()
		{ return _col; }

	bool operator!= (const ColIterator& c) const
		{ return (_col.begin () != c._col.begin ()) || (_col.end () != c._col.end ()); }
   
	operator ConstColIterator ()
	{
            return ConstColIterator (reinterpret_cast<ConstCol&> (_col) , _stride);
	}

    private:

	Col _col;
	size_t _stride;
};

template <class _Element>
template <class Field>
DenseMatrixBase<_Element>::DenseMatrixBase( MatrixStream<Field>& ms )
	:_rep(0), _rows(0), _cols(0), _ptr(NULL)
{
	if( !ms.getArray(_rep) || !ms.getRows(_rows) || !ms.getColumns(_cols) )
		throw ms.reportError(__FUNCTION__,__LINE__);
	_ptr = &_rep[0];
}

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
typename DenseMatrixBase<Element>::RowIterator DenseMatrixBase<Element>::rowBegin ()
	{ return RowIterator (_rep.begin (), _cols, _cols); }

template <class Element>
typename DenseMatrixBase<Element>::RowIterator DenseMatrixBase<Element>::rowEnd ()
	{ return RowIterator (_rep.end (), _cols, _cols); }
  
template <class Element>
typename DenseMatrixBase<Element>::ConstRowIterator DenseMatrixBase<Element>::rowBegin () const
	{ return ConstRowIterator (_rep.begin (), _cols, _cols); }  

template <class Element>
typename DenseMatrixBase<Element>::ConstRowIterator DenseMatrixBase<Element>::rowEnd () const
	{return ConstRowIterator (_rep.end (), _cols, _cols); }
  
template <class Element>
typename DenseMatrixBase<Element>::ColIterator DenseMatrixBase<Element>::colBegin ()
	{ return  typename DenseMatrixBase<Element>::ColIterator (_rep.begin (), _cols, _rows); }

template <class Element>
typename DenseMatrixBase<Element>::ColIterator DenseMatrixBase<Element>::colEnd ()
	{ return  typename DenseMatrixBase<Element>::ColIterator (_rep.begin ()+_cols, _cols, _rows); }
  
template <class Element>
typename DenseMatrixBase<Element>::ConstColIterator DenseMatrixBase<Element>::colBegin () const
	{ return  typename DenseMatrixBase<Element>::ConstColIterator (_rep.begin (), _cols, _rows); }
  
template <class Element>
typename DenseMatrixBase<Element>::ConstColIterator DenseMatrixBase<Element>::colEnd () const
	{ return  typename DenseMatrixBase<Element>::ConstColIterator (_rep.begin ()+_cols, _cols, _rows); }

template <class Element>
template <class Field>
std::istream &DenseMatrixBase<Element>::read (std::istream &file, const Field &F)
{
//  	RawIterator p;
//  	int m,n; 
//  	char c;
//  	file>>m>>n>>c;
	
//  	if (m*n < _rows*_cols)
//  		cerr<<"NOT ENOUGH ELEMENT TO READ\n";
// 	else {
//  		for (p = rawBegin (); p != rawEnd (); ++p) {	
// 			integer tmp;
// 			file>>tmp;cout<<tmp<<endl;
//  			//file.ignore(1); 	
// 			F.read (file, *p);
//  		}
//  	}


  	RawIterator p;
  	int m,n; 
  	char c;
	file>>m>>n>>c;
  	_rows = m; _cols = n; 

	Element zero;
	F.init(zero,0UL);
	_rep.resize(_rows * _cols, zero);
	_ptr= &_rep[0];
	 	
	if ((c != 'M') && (c != 'm'))
 	 	for (p = rawBegin (); p != rawEnd (); ++p) 
			{	
				//file.ignore(1);
				F.read (file, *p);
			}

	else { // sparse file format - needs fixing
		int i, j;
		while (true)
			{	
				file >> i >> j;
				//file.ignore(1);
				//if (! file) break;
				if (i+j <= 0) break;
				F.read (file, _rep[_cols*(i-1) + j-1]);			
			}
	}
	
	return file;	
}
  
template <class Element>
template <class Field>
std::ostream& DenseMatrixBase<Element>::write (std::ostream &os, const Field &F) const
{

	integer c;
	int wid;

	ConstRowIterator p;



	F.cardinality (c);

	if (c >0)
		wid = (int) ceil (log ((double) c) / M_LN10);
	else {
		integer tmp;
		size_t max=0;
		ConstRawIterator it = rawBegin();
		for (; it != rawEnd(); ++it){
			F.convert(tmp,*it);
			if (tmp.bitsize() > max)
				max= tmp.bitsize();
		}
		wid= (int) ceil ((double)max / M_LN10)+1; 
	}

	for (p = rowBegin (); p != rowEnd (); ++p) {
		typename ConstRow::const_iterator pe;

		os << "  [ ";

		for (pe = p->begin (); pe != p->end (); ++pe) {
			os.width (wid);
			F.write (os, *pe);
			os << " ";
		}

		os << "]" << std::endl;
	}

	return os;
}

template <class Element>
class DenseMatrixBase<Element>::RawIndexedIterator
{
        size_t _r_index;
	size_t _c_index;
	size_t _dim;
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

	
	RawIndexedIterator operator ++ (int)
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


	RawIndexedIterator operator -- (int) 
	{
		RawIndexedIterator tmp = *this;
		--(*this);
		return tmp;
	}	

	Element &operator * () const
		{ return *(_begin + (_r_index * _dim + _c_index)); }

	
	Element* operator -> () const
		{ return _begin + (_r_index * _dim + _c_index); }
	

	size_t rowIndex () const
		{ return _r_index; }

	size_t colIndex () const
	{ return _c_index; }

	const Element &value () const
        { return *(_begin + (_r_index * _dim + _c_index)); }

	
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
class DenseMatrixBase<Element>::ConstRawIndexedIterator
{
	size_t _r_index;
	size_t _c_index;
	size_t _dim;
    	typedef Element value_type;
	typename Rep::const_iterator _begin;

    public:
	ConstRawIndexedIterator (const size_t  &dim,
			    const size_t  &r_index,
			    const size_t  &c_index,
			    const typename Rep::const_iterator &begin)
		: _r_index (r_index), _c_index (c_index), _dim (dim), _begin (begin)
	{}
	
	ConstRawIndexedIterator ():_r_index (0), _c_index (0), _dim (1), _begin (0){}

	ConstRawIndexedIterator (const ConstRawIndexedIterator& r)
		: _r_index (r._r_index), _c_index (r._c_index), _dim (r._dim), _begin (r._begin)
	{}

	ConstRawIndexedIterator& operator = (const ConstRawIndexedIterator &iter)
	{
		_r_index = iter._r_index;
		_c_index = iter._c_index;
		_dim = iter._dim;
		_begin = iter._begin;
		return *this;
	}
	
	bool operator == (const ConstRawIndexedIterator &iter) const
		{ return (_r_index == iter._r_index) &&
			  (_c_index == iter._c_index) &&
			  (_dim == iter._dim) &&
			  (_begin==iter._begin); }

	bool operator != (const ConstRawIndexedIterator& iter) const
		{ return (_r_index != iter._r_index) ||
			  (_c_index != iter._c_index) ||
			  (_dim != iter._dim) ||
			  (_begin!=iter._begin); }
	
	ConstRawIndexedIterator &operator ++ ()
	{
		++_c_index;

		if (_c_index == _dim) {
			_c_index = 0;
			++_r_index;
		}

		return *this;
	}

	
	ConstRawIndexedIterator operator ++ (int) 
	{
		ConstRawIndexedIterator tmp = *this;
		++(*this);
		return tmp;
	}

	ConstRawIndexedIterator &operator -- ()
	{ 
		if (_c_index)
			--_c_index;
		else {
			--_r_index;
			_c_index = _dim - 1;
		}

		return *this;
	}



	ConstRawIndexedIterator operator -- (int) 
	{
		ConstRawIndexedIterator tmp = *this;
		--(*this);
		return tmp;
	}	

	const Element &operator * () const
		{ return *(_begin + (_r_index * _dim + _c_index)); }

	const Element *operator -> () const
		{ return _begin + (_r_index * _dim + _c_index); }

	size_t rowIndex () const
		{ return _r_index; }

	size_t colIndex () const
		{ return _c_index; }
	
    	const Element &value() const
		{ return *(_begin + (_r_index * _dim + _c_index)); }
};



template <class Element>
typename DenseMatrixBase<Element>::ConstRawIndexedIterator DenseMatrixBase<Element>::rawIndexedBegin () const
{
	return ConstRawIndexedIterator (coldim (), 0, 0, _rep.begin ());
}

template <class Element>
typename DenseMatrixBase<Element>::ConstRawIndexedIterator DenseMatrixBase<Element>::rawIndexedEnd () const
{
	return ConstRawIndexedIterator (coldim (), rowdim (), 0, _rep.begin ());
}

} // namespace LinBox

#endif // __LINBOX_matrix_dense_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
