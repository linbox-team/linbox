/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/sparse.inl
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
 *   - Eliminated SparseMatrixAux and moved that functionality into Sparse0
 *   - Made SparseMatrixBase parameterized only on the element type
 *   - New read/write implementations for SparseMatrixBase, supporting multiple
 *     formats
 *   - Eliminated Gaussian elimination code
 *   - Added iterators, including ColOfRowsIterator, RawIterator, and
 *     RawIndexIterator
 *   - Eliminated operator []; added getEntry; changed put_value to setEntry
 * ------------------------------------
 * 
 * See COPYING for license information.
 */

#ifndef __BLACKBOX_SPARSE_INL
#define __BLACKBOX_SPARSE_INL

#include "linbox/blackbox/sparse.h"

namespace LinBox 
{

template <class Field, class _Vector, class _Row, class VectorTrait>
_Vector &SparseMatrix<Field, _Vector, _Row, VectorCategories::DenseVectorTag<VectorTrait> >::apply
	(_Vector &y, const _Vector &x) const
{
	linbox_check (x.size () == _n);

	typename std::vector<Row>::const_iterator i;
	typename Vector::iterator y_iter = y.begin();
		
	for (i = _A.begin (); i != _A.end (); i++, y_iter++)
		_VD.dot (*y_iter, *i, x);

	return y;
}

template <class Field, class _Vector, class _Row, class VectorTrait>
inline _Vector &SparseMatrix<Field, _Vector, _Row, VectorCategories::DenseVectorTag<VectorTrait> >::applyTransposeSpecialized
	(_Vector &y, const _Vector &x) const
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

template <class Field, class _Vector, class _Row, class VectorTrait>
template <class RowTrait>
inline _Vector &SparseMatrix<Field, _Vector, _Row, VectorCategories::DenseVectorTag<VectorTrait> >::applyTransposeSpecialized
	(_Vector &y, const _Vector &x, VectorCategories::SparseParallelVectorTag<RowTrait> tag) const
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
		typename Row::first_type::const_iterator k_idx;
		typename Row::second_type::const_iterator k_elt;
		typename Vector::const_iterator j;

		for (i = _A.begin (), j = x.begin (); i != _A.end (); i++, j++)
			for (k_idx = i->first.begin (), k_elt = i->second.begin (); k_idx != i->first.end (); ++k_idx, ++k_elt)
				_faxpy[*k_idx].accumulate (*k_elt, *j);
	}

	{
		typename Vector::iterator i;
		typename std::vector<FieldAXPY<Field> >::iterator j;

		for (i = y.begin (), j = _faxpy.begin (); j != _faxpy.end (); i++, j++)
			j->get (*i);
	}

	return y;
}
 
template <class Field, class _Vector, class _Row, class VectorTrait>
_Vector &SparseMatrix<Field, _Vector, _Row, VectorCategories::SparseSequenceVectorTag<VectorTrait> >::apply
	(_Vector &y, const _Vector &x) const
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

template <class Field, class _Vector, class _Row, class VectorTrait>
inline _Vector &SparseMatrix<Field, _Vector, _Row, VectorCategories::SparseSequenceVectorTag<VectorTrait> >::applyTransposeSpecialized
	(_Vector &y, const _Vector &x) const
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

template <class Field, class _Vector, class _Row, class VectorTrait>
template <class RowTrait>
inline _Vector &SparseMatrix<Field, _Vector, _Row, VectorCategories::SparseSequenceVectorTag<VectorTrait> >::applyTransposeSpecialized
	(_Vector &y, const _Vector &x, VectorCategories::SparseParallelVectorTag<RowTrait> tag) const
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
		typename Row::first_type::const_iterator k_idx;
		typename Row::second_type::const_iterator k_elt;
		typename Vector::const_iterator j;

		for (j = x.begin (); j != x.end (); j++)
			for (k_idx = _A[j->first].first.begin (), k_elt = _A[j->first].second.begin ();
			     k_idx != _A[j->first].first.end ();
			     ++k_idx, ++k_elt)
				_faxpy[*k_idx].accumulate (*k_elt, j->second);
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

template <class Field, class _Vector, class _Row, class VectorTrait>
_Vector &SparseMatrix<Field, _Vector, _Row, VectorCategories::SparseAssociativeVectorTag<VectorTrait> >::apply
	(_Vector &y, const _Vector &x) const
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

template <class Field, class _Vector, class _Row, class VectorTrait>
inline _Vector &SparseMatrix<Field, _Vector, _Row, VectorCategories::SparseAssociativeVectorTag<VectorTrait> >::applyTransposeSpecialized
	(_Vector &y, const _Vector &x) const
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

template <class Field, class _Vector, class _Row, class VectorTrait>
template <class RowTrait>
inline _Vector &SparseMatrix<Field, _Vector, _Row, VectorCategories::SparseAssociativeVectorTag<VectorTrait> >::applyTransposeSpecialized
	(_Vector &y, const _Vector &x, VectorCategories::SparseParallelVectorTag<RowTrait> tag) const
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
		typename Row::first_type::const_iterator k_idx;
		typename Row::second_type::const_iterator k_elt;
		typename Vector::const_iterator j;

		for (j = x.begin (); j != x.end (); j++)
			for (k_idx = _A[j->first].first.begin (), k_elt = _A[j->first].second.begin ();
			     k_idx != _A[j->first].first.end ();
			     ++k_idx, ++k_elt)
				_faxpy[*k_idx].accumulate (*k_elt, j->second);
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
 
template <class Field, class _Vector, class _Row, class VectorTrait>
_Vector &SparseMatrix<Field, _Vector, _Row, VectorCategories::SparseParallelVectorTag<VectorTrait> >::apply
	(_Vector &y, const _Vector &x) const
{
	linbox_check ((x.first.size () == 0) || (x.first.back () < _n));

	y.first.clear ();
	y.second.clear ();

	typename std::vector<Row>::const_iterator i;
	int idx;
	Element tmp;

	for (i = _A.begin (), idx = 0; i != _A.end (); i++, idx++) {
		_VD.dot (tmp, *i, x);

		if (!_F.isZero (tmp)) {
			y.first.push_back (idx);
			y.second.push_back (tmp);
		}
	}
    
	return y;
}

template <class Field, class _Vector, class _Row, class VectorTrait>
inline _Vector &SparseMatrix<Field, _Vector, _Row, VectorCategories::SparseParallelVectorTag<VectorTrait> >::applyTransposeSpecialized
	(_Vector &y, const _Vector &x) const
{
	y.first.clear ();
	y.second.clear ();

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
		typename Vector::first_type::const_iterator j_idx;
		typename Vector::second_type::const_iterator j_elt;

		for (j_idx = x.first.begin (), j_elt = x.second.begin (); j_idx != x.first.end (); ++j_idx, ++j_elt)
			for (k = _A[*j_idx].begin (); k != _A[*j_idx].end (); ++k)
				_faxpy[k->first].accumulate (k->second, *j_elt);
	}
		
	{
		size_t i;
		Element tmp;
		typename std::vector<FieldAXPY<Field> >::iterator j;

		for (j = _faxpy.begin (), i = 0; j != _faxpy.end (); j++, i++) {
			j->get (tmp);
			if (!_F.isZero (tmp)) {
				y.first.push_back (i);
				y.second.push_back (tmp);
			}
		}
	}

	return y;
}

template <class Field, class _Vector, class _Row, class VectorTrait>
template <class RowTrait>
inline _Vector &SparseMatrix<Field, _Vector, _Row, VectorCategories::SparseParallelVectorTag<VectorTrait> >::applyTransposeSpecialized
	(_Vector &y, const _Vector &x, VectorCategories::SparseParallelVectorTag<RowTrait> tag) const
{
	y.first.clear ();
	y.second.clear ();

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
		typename Row::first_type::const_iterator k_idx;
		typename Row::second_type::const_iterator k_elt;
		typename Vector::first_type::const_iterator j_idx;
		typename Vector::second_type::const_iterator j_elt;

		for (j_idx = x.first.begin (), j_elt = x.second.begin (); j_idx != x.first.end (); ++j_idx, ++j_elt)
			for (k_idx = _A[*j_idx].first.begin (), k_elt = _A[*j_idx].second.begin ();
			     k_idx != _A[*j_idx].first.end ();
			     ++k_idx, ++k_elt)
				_faxpy[*k_idx].accumulate (*k_elt, *j_elt);
	}
		
	{
		size_t i;
		Element tmp;
		typename std::vector<FieldAXPY<Field> >::iterator j;

		for (j = _faxpy.begin (), i = 0; j != _faxpy.end (); j++, i++) {
			j->get (tmp);

			if (!_F.isZero (tmp)) {
				y.first.push_back (i);
				y.second.push_back (tmp);
			}
		}
	}

	return y;
}

template <class Field, class BElement, class _Vector, class _Row, class BRow>
BlackboxArchetype<_Vector> *SparseMatrixFactory<Field, BElement, _Vector, _Row, BRow>::makeBlackbox (const Field &F)
{
	SparseMatrix<Field, _Vector, _Row> *A = new SparseMatrix<Field, _Vector, _Row> (F, rowdim (), coldim ());

	typename SparseMatrixBase<BElement, BRow>::ConstRawIterator i;
	typename SparseMatrixBase<BElement, BRow>::ConstRawIndexedIterator j;

	for (i = _A.rawBegin (), j = _A.rawIndexedBegin (); i != _A.rawEnd (); ++i, ++j)
		F.init (A->refEntry (j.rowIndex (), j.colIndex ()), *i);

	return A;
}

}

#endif // __BLACKBOX_SPARSE_INL
