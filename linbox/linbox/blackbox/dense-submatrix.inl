/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/dense-submatrix.inl
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

#ifndef DENSE_SUBMATRIX_INL
#define DENSE_SUBMATRIX_INL

#include "linbox/util/debug.h"
#include "linbox/blackbox/dense.h"
#include "linbox/blackbox/dense-submatrix.h"

namespace LinBox
{
 
template <class Field, class Vector>
DenseSubmatrix<Field, Vector>::DenseSubmatrix (DenseMatrix<Field, Vector> *M,
					       size_t row,
					       size_t col,
					       size_t rowdim,
					       size_t coldim)
	: _M (M), _beg_row (row), _end_row (row + rowdim), _beg_col (col), _end_col (col + coldim)
{
	linbox_check (_beg_row < _end_row);
	linbox_check (_beg_col < _end_col);
	linbox_check (_end_row <= M->rowdim ());
	linbox_check (_end_col <= M->coldim ());
} 

template <class Field, class Vector>
DenseSubmatrix<Field, Vector>::DenseSubmatrix (const DenseSubmatrix<Field, Vector> &SM,
					       size_t row,
					       size_t col,
					       size_t rowdim,
					       size_t coldim)
	: _M (SM._M),
	  _beg_row (SM._beg_row + row),
	  _end_row (SM._beg_row + rowdim),
	  _beg_col (SM._beg_col + col),
	  _end_col (SM._beg_col + coldim)
{
	linbox_check (_beg_row < _end_row);
	linbox_check (_beg_col < _end_col);
	linbox_check (_end_row - _beg_row < SM.rowdim ());
	linbox_check (_end_col - _beg_col < SM.coldim ());
}
  
template <class Field, class Vector>
DenseSubmatrix<Field, Vector>::DenseSubmatrix (const DenseSubmatrix<Field, Vector> &SM)
	: _M (SM._M), _beg_row (SM._beg_row), _end_row (SM._end_row), _beg_col (SM._beg_col), _end_col (SM._end_col)
{}

template <class Field, class Vector>
DenseSubmatrix<Field, Vector>& DenseSubmatrix<Field, Vector>::operator=(const DenseSubmatrix<Field, Vector> &SM)
{
	_M = SM._M;
	_beg_row = SM._beg_row;
	_end_row = SM._end_row;
	_beg_col = SM._beg_col;
	_end_col = SM._end_col;

	return *this;
}
  
template <class Field, class Vector>
class DenseSubmatrix<Field, Vector>::RawIterator
{
    public:
	RawIterator (){}

	RawIterator (const typename DenseMatrix<Field>::RawIterator& beg,
		     const typename DenseMatrix<Field>::RawIterator& cur,
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
		RawIterator tmp = *this;
		this->operator++ ();
		return tmp;
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
	typename DenseMatrix<Field>::RawIterator _cur;
	typename DenseMatrix<Field>::RawIterator _beg;
	size_t _cont_len;
	size_t _gap_len;
};
  
template <class Field, class Vector>
class DenseSubmatrix<Field, Vector>::ConstRawIterator
{   
    public:
	ConstRawIterator (){}

	ConstRawIterator (const typename DenseMatrix<Field>::ConstRawIterator& beg, 
			  const typename DenseMatrix<Field>::ConstRawIterator& cur, 
			  size_t cont_len,
			  size_t gap_len)
		: _beg (beg), _cur (cur), _cont_len (cont_len), _gap_len (gap_len)
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

	ConstRawIterator& operator++(int)
	{
		RawIterator tmp = *this;
		this->operator++();
		return tmp;
	}

	bool operator != (const ConstRawIterator& r) const
	{
		return (_cur != r._cur) || (_beg != r._beg) || (_cont_len != r._cont_len) || (_gap_len != r._gap_len);
	}
    
	Element& operator*()
		{ return *_cur; }

	const Element& operator*() const
		{ return *_cur; }

    protected:
	typename DenseMatrix<Field>::ConstRawIterator _cur;
	typename DenseMatrix<Field>::ConstRawIterator _beg;
	size_t _cont_len;
	size_t _gap_len;
};
  
template <class Field, class Vector>
typename DenseSubmatrix<Field, Vector>::RawIterator DenseSubmatrix<Field, Vector>::rawBegin ()
{
	return RawIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col, 
			    _M->rawBegin () + _beg_row * _M->coldim () + _beg_col, 
			    coldim (), _M->coldim () - coldim ());
}
  
template <class Field, class Vector>
typename DenseSubmatrix<Field, Vector>::RawIterator DenseSubmatrix<Field, Vector>::rawEnd ()
{
	return RawIterator (_M->rawBegin () + _end_row * _M->coldim () + _beg_col, 
			    _M->rawBegin () + _end_row * _M->coldim () + _beg_col, 
			    coldim (), _M->coldim () - coldim ());
}
    
template <class Field, class Vector>
typename DenseSubmatrix<Field, Vector>::ConstRawIterator DenseSubmatrix<Field, Vector>::rawBegin () const
{
	return ConstRawIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col, 
				 _M->rawBegin () + _beg_row * _M->coldim () + _beg_col, 
				 coldim (), _M->coldim () - coldim ());
}
  
template <class Field, class Vector>
typename DenseSubmatrix<Field, Vector>::ConstRawIterator DenseSubmatrix<Field, Vector>::rawEnd () const
{
	return ConstRawIterator (_M->rawBegin () + _end_row * _M->coldim () + _beg_col, 
				 _M->rawBegin () + _end_row * _M->coldim () + _beg_col, 
				 coldim (), _M->coldim () - coldim ());
}
  
template <class Field, class Vector>
typename DenseSubmatrix<Field, Vector>::ColOfRowsIterator DenseSubmatrix<Field, Vector>::colOfRowsBegin ()
{
	return ColOfRowsIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col,
				  _end_col - _beg_col, _M->coldim ());
}
 
template <class Field, class Vector>
typename DenseSubmatrix<Field, Vector>::ColOfRowsIterator DenseSubmatrix<Field, Vector>::colOfRowsEnd ()
{
	return ColOfRowsIterator (_M->rawBegin () + _end_row * _M->coldim () + _beg_col,
				  _end_col - _beg_col, _M->coldim ());
}

template <class Field, class Vector>
typename DenseSubmatrix<Field, Vector>::ConstColOfRowsIterator DenseSubmatrix<Field, Vector>::colOfRowsBegin () const
{
	return ConstColOfRowsIterator (_M->rawBegin () + _beg_row * _M->coldim () + _beg_col,
				       _end_col - _beg_col, _M->coldim ());
}
  
template <class Field, class Vector>
typename DenseSubmatrix<Field, Vector>::ConstColOfRowsIterator DenseSubmatrix<Field, Vector>::colOfRowsEnd () const
{
	return ConstColOfRowsIterator (_M->rawBegin () + _end_row * _M->coldim () + _beg_col,
				       _end_col - _beg_col, _M->coldim ());
}

template <class Field, class Vector>
typename DenseSubmatrix<Field, Vector>::RowOfColsIterator DenseSubmatrix<Field, Vector>::rowOfColsBegin ()
{
	return RowOfColsIterator (_M->rawBegin () + _beg_col + _beg_row * _M->coldim (),
				  _M->coldim (), rowdim ());
}

template <class Field, class Vector>
typename DenseSubmatrix<Field, Vector>::RowOfColsIterator DenseSubmatrix<Field, Vector>::rowOfColsEnd ()
{
	return RowOfColsIterator (_M->rawBegin () + _end_col + _beg_row * _M->coldim (),
				  _M->coldim (), rowdim ());
}

template <class Field, class Vector>
typename DenseSubmatrix<Field, Vector>::ConstRowOfColsIterator DenseSubmatrix<Field, Vector>::rowOfColsBegin () const
{
	return ConstRowOfColsIterator (_M->rawBegin () + _beg_col + _beg_row * _M->coldim (),
				       _M->coldim (), rowdim ());
}

template <class Field, class Vector>
typename DenseSubmatrix<Field, Vector>::ConstRowOfColsIterator DenseSubmatrix<Field, Vector>::rowOfColsEnd () const
{
	return ConstRowOfColsIterator (_M->rawBegin () + _end_col + _beg_row * _M->coldim (),
				       _M->coldim (), rowdim ());
}

template <class Field, class Vector>
void DenseSubmatrix<Field, Vector>::read (std::istream &file)
{
	RawIterator p;

	for (p = rawBegin (); p != rawEnd (); ++p) {
		file.ignore (1);
		_M->field ().read (file, *p);
	}
}

template <class Field, class Vector>
std::ostream &DenseSubmatrix<Field, Vector>::write (std::ostream &os) const
{
	ConstColOfRowsIterator p;

	for (p = colOfRowsBegin (); p != colOfRowsEnd (); ++p) {
		ConstRowIterator pe;

		for (pe = p->begin (); pe != p->end (); ++pe) {
			_M->field ().write (os, *pe);
			os << ' ';
		}

		os << endl;
	}

	os << endl;

	return os;
}

template <class Field, class Vector>
template<class Vect1, class Vect2>
Vect1& DenseSubmatrix<Field, Vector>::apply (Vect1& y, const Vect2& x) const
{
	ConstColOfRowsIterator         p;
	typename Vect1::iterator       p_y = y.begin ();  
	typename Vect2::const_iterator p_x;

	for (p = colOfRowsBegin (); p != colOfRowsEnd (); ++p, ++p_y) {
		_M->field ().init (*p_y, 0);

		for (ConstRowIterator pe = p->begin (), p_x = x.begin (); pe != p->end (); ++pe, ++p_x)
			_M->field ().axpyin (*p_y, *pe, *p_x);
	}
    
	return y;
}
  
template <class Field, class Vector>
template<class Iterator1, class Iterator2>
Iterator1& DenseSubmatrix<Field, Vector>::apply (Iterator1 in,
						 const Iterator2& outbegin,
						 const Iterator2& outend) const
{
	linbox_check (coldim () == (outend - outbegin));

	ConstColOfRowsIterator rowp;
	Iterator2 p_out;
	ConstRowIterator pe;

	for (rowp = colOfRowsBegin (); rowp != colOfRowsEnd (); ++rowp, ++in) {
		_M->field ().init (*in, 0);

		for (pe = rowp->begin (), p_out = outbegin; pe != rowp->end (); ++pe, ++p_out)
			_M->field ().axpyin (*in, *pe, *p_out);
	}
    
	return in;
}

template <class Field, class Vector>
template <class Vect1, class Vect2>
Vect1& DenseSubmatrix<Field, Vector>::applyTranspose (Vect1& y, const Vect2& x) const
{
	ConstRowOfColsIterator colp;
	typename Vect1::iterator  p_y = y.begin ();  
	typename Vect2::const_iterator p_x;

	for (colp = rowOfColsBegin (); colp != rowOfColsEnd (); ++colp, ++p_y) {
		_M->field ().init (*p_y, 0);

		for (ConstColIterator pe = colp->begin (), p_x = x.begin (); pe != colp->end (); ++pe, ++p_x)
			_M->field ().axpyin (*p_y, *pe, *p_x);
	}
    
	return y;
}

template <class Field, class Vector>
template <class Iterator1, class Iterator2>
Iterator1& DenseSubmatrix<Field, Vector>::applyTranspose (Iterator1 in,
							  const Iterator2& outbegin,
							  const Iterator2& outend) const
{
	linbox_check (rowdim () == (outend - outbegin));

	ConstRowOfColsIterator colp;
	Iterator2 p_out;
	ConstColIterator pe;

	for (colp = rowOfColsBegin (); colp != rowOfColsEnd (); ++colp, ++in) {
		_M->field ().init (*in, 0);

		for (pe = colp->begin (), p_out = outbegin; pe != colp->end (); ++pe, ++p_out)
			_M->field ().axpyin (*in, *pe, *p_out);
	}
    
	return in;
}

template <class Field, class Vector>
typename DenseSubmatrix<Field, Vector>::Row DenseSubmatrix<Field, Vector>::operator [] (int i)
	{ return Row (_M->rawBegin () + i * _M->coldim () + _beg_row, 
		      _M->rawBegin () + i * _M->coldim () + _end_row); }

template <class Field, class Vector>
typename DenseSubmatrix<Field, Vector>::ConstRow DenseSubmatrix<Field, Vector>::operator[](int i) const
	{ return ConstRow (_M->rawBegin () + i * _M->coldim () + _beg_row, 
			   _M->rawBegin () + i * _M->coldim () + _end_row); }
    

} // namespace LinBox

#endif // DENSE_SUBMATRIX_INL
