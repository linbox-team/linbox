/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/sparse0.C
 * Copyright (C) 1999-2001 William J Turner,
 *               2001-2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * Modified by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Refactoring:
 *   - Eliminated SparseMatrix0Aux and moved that functionality into Sparse0
 *   - Made SparseMatrix0Base parameterized only on the element type
 *   - New read/write implementations for SparseMatrix0Base, supporting multiple
 *     formats
 *   - Eliminated Gaussian elimination code
 *   - Added iterators, including ColOfRowsIterator, RawIterator, and
 *     RawIndexIterator
 *   - Eliminated operator []; added getEntry; changed put_value to setEntry
 * ------------------------------------
 * 
 * See COPYING for license information.
 */

#ifndef __SPARSE0_C
#define __SPARSE0_C

#include "linbox/blackbox/sparse0.h"

namespace LinBox 
{

template <class Field, class Row, class Vector, class VectorTrait>
Vector &SparseMatrix0<Field, Vector, Row, VectorCategories::DenseVectorTag<VectorTrait> >::apply
	(Vector &y, const Vector &x) const
{
	linbox_check (x.size () == _n);

	typename std::vector<Row>::const_iterator i;
	typename Vector::iterator y_iter = y.begin();
		
	for (i = _A.begin (); i != _A.end (); i++, y_iter++)
		_VD.dot (*y_iter, *i, x);

	return y;
}

template <class Field, class Row, class Vector, class VectorTrait>
Vector &SparseMatrix0<Field, Vector, Row, VectorCategories::DenseVectorTag<VectorTrait> >::applyTranspose
	(Vector &y, const Vector &x) const
{
	linbox_check (x.size () == _m);

	if (_faxpy.size () == 0) {
		for (int i = _n; i--;)
			_faxpy.push_back (FieldAXPY<Field> (_F));
	} else {
		typename Field::Element zero;
		typename std::vector<FieldAXPY<Field> >::iterator i;

		_F.init (zero, 0);

		for (i = _faxpy.begin (); i != _faxpy.end (); i++)
			i->assign (zero);
	}

	{
		typename Rep::const_iterator i;
		typename Row::const_iterator k;
		typename Vector::const_iterator j;

		for (i = _A.begin (), j = x.begin (); i != _A.end (); i++, j++)
			for (k = i->begin (); k != i->end (); k++)
				_faxpy[k->first].accumulate (k->second, *j);
	}
		
	{
		typename Vector::iterator i;
		typename std::vector<FieldAXPY<Field> >::iterator j;

		for (i = y.begin (), j = _faxpy.begin (); j != _faxpy.end (); i++, j++)
			j->get (*i);
	}

	return y;
}
 
template <class Field, class Row, class Vector, class VectorTrait>
Vector &SparseMatrix0<Field, Vector, Row, VectorCategories::SparseSequenceVectorTag<VectorTrait> >::apply
	(Vector &y, const Vector &x) const
{
	linbox_check ((x.size () == 0) || (x.back ().first < _n));

	y.clear ();

	typename std::vector<Row>::const_iterator i;
	int idx;
	Element tmp;

	for (i = _A.begin (), idx = 0; i != _A.end (); i++, idx++) {
		_VD.dot (tmp, *i, x);
		if (!_F.isZero (tmp))
			y.push_back (std::pair <size_t, typename Field::Element> (idx, tmp));
	}
    
	return y;
}

template <class Field, class Row, class Vector, class VectorTrait>
Vector &SparseMatrix0<Field, Vector, Row, VectorCategories::SparseSequenceVectorTag<VectorTrait> >::applyTranspose
	(Vector &y, const Vector &x) const
{
	y.clear ();

	if (_faxpy.size () == 0) {
		for (int i = _n; i--;)
			_faxpy.push_back (FieldAXPY<Field> (_F));
	} else {
		typename Field::Element zero;
		typename std::vector<FieldAXPY<Field> >::iterator i;

		_F.init (zero, 0);

		for (i = _faxpy.begin (); i != _faxpy.end (); i++)
			i->assign (zero);
	}

	{
		typename Row::const_iterator k;
		typename Vector::const_iterator j;

		for (j = x.begin (); j != x.end (); j++)
			for (k = _A[j->first].begin (); k != _A[j->first].end (); k++)
				_faxpy[k->first].accumulate (k->second, j->second);
	}
		
	{
		size_t i;
		Element tmp;
		typename std::vector<FieldAXPY<Field> >::iterator j;

		for (j = _faxpy.begin (), i = 0; j != _faxpy.end (); j++, i++) {
			j->get (tmp);
			if (!_F.isZero (tmp))
				y.push_back (std::pair<size_t, Element> (i, tmp));
		}
	}

	return y;
}

template <class Field, class Row, class Vector, class VectorTrait>
Vector &SparseMatrix0<Field, Vector, Row, VectorCategories::SparseAssociativeVectorTag<VectorTrait> >::apply
	(Vector &y, const Vector &x) const
{
	linbox_check ((x.size () == 0) || (x.rbegin ()->first < _n));
 
	y.clear ();

	typename std::vector<Row>::const_iterator i;
	size_t idx;
	Element tmp;
 
	for (i = _A.begin (), idx = 0; i != _A.end (); i++, idx++) {
		_VD.dot (tmp, *i, x);
		if (!_F.isZero (tmp))
			y[idx] = tmp;
	}
    
	return y;
}

template <class Field, class Row, class Vector, class VectorTrait>
Vector &SparseMatrix0<Field, Vector, Row, VectorCategories::SparseAssociativeVectorTag<VectorTrait> >::applyTranspose
	(Vector &y, const Vector &x) const
{
	y.clear ();

	if (_faxpy.size () == 0) {
		for (int i = _n; i--;)
			_faxpy.push_back (FieldAXPY<Field> (_F));
	} else {
		typename Field::Element zero;
		typename std::vector<FieldAXPY<Field> >::iterator i;

		_F.init (zero, 0);

		for (i = _faxpy.begin (); i != _faxpy.end (); i++)
			i->assign (zero);
	}

	{
		typename Row::const_iterator k;
		typename Vector::const_iterator j;

		for (j = x.begin (); j != x.end (); j++)
			for (k = _A[j->first].begin (); k != _A[j->first].end (); k++)
				_faxpy[k->first].accumulate (k->second, j->second);
	}
		
	{
		size_t i;
		Element tmp;
		typename std::vector<FieldAXPY<Field> >::iterator j;

		for (j = _faxpy.begin (), i = 0; j != _faxpy.end (); j++, i++) {
			j->get (tmp);
			if (!_F.isZero (tmp))
				y[i] = tmp;
		}
	}

	return y;
}
 
}

#endif // __SPARSE0_C
