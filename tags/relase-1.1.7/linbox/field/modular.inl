/* linbox/field/modular.inl
 * Copyright (C) 2002 Bradford Hovinen
 * Copyright (C) 2002 Ahmet Duran
 * Copyright (C) 2002 B. David Saunders
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 *            Ahmet Duran <duran@cis.udel.edu>,
 *            Dave Saunders <saunders@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_field_modular_INL
#define __LINBOX_field_modular_INL

//Dan Roche 7-2-04
#ifndef __LINBOX_MIN
#define __LINBOX_MIN(a,b) ( (a) < (b) ? (a) : (b) )
#endif

#include <iostream>

namespace LinBox {

template <class Vector1, class Vector2>
inline uint8 &DotProductDomain<Modular<uint8> >::dotSpecializedDD
	(uint8 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::const_iterator i = v1.begin ();
	typename Vector2::const_iterator j = v2.begin ();

	typename Vector1::const_iterator iterend = v1.begin () + v1.size() % _F._k;

	uint64 y = 0;

	for (; i != iterend; ++i, ++j)
		y += (uint64) *i * (uint64) *j;

	y %= (uint64) _F._modulus;

	for (; iterend != v1.end (); j += _F._k) {
		typename Vector1::const_iterator iter_i = iterend;
		typename Vector2::const_iterator iter_j;

		iterend += _F._k;

		for (iter_j = j; iter_i != iterend; ++iter_i, ++iter_j)
			y += (uint64) *iter_i * (uint64) *j;

		y %= (uint64) _F._modulus;
	}

	return res = y;
}

template <class Vector1, class Vector2>
inline uint8 &DotProductDomain<Modular<uint8> >::dotSpecializedDSP
	(uint8 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::first_type::const_iterator i_idx = v1.first.begin ();
	typename Vector1::second_type::const_iterator i_elt = v1.second.begin ();

	uint64 y = 0;

	if (v1.first.size () < _F._k) {
		for (; i_idx != v1.first.end (); ++i_idx, ++i_elt)
			y += (uint64) *i_elt * (uint64) v2[*i_idx];

		return res = y % (uint64) _F._modulus;
	} else {
		typename Vector1::first_type::const_iterator iterend = v1.first.begin () + v1.first.size() % _F._k;

		for (; i_idx != iterend; ++i_idx, ++i_elt)
			y += (uint64) *i_elt * (uint64) v2[*i_idx];

		y %= (uint64) _F._modulus;

		while (iterend != v1.first.end ()) {
			typename Vector1::first_type::const_iterator iter_i_idx = iterend;
			typename Vector1::second_type::const_iterator iter_i_elt = i_elt;

			iterend += _F._k;
			i_elt += _F._k;

			for (; iter_i_idx != iterend; ++iter_i_idx, ++iter_i_elt)
				y += (uint64) *iter_i_elt * (uint64) v2[*iter_i_idx];

			y %= (uint64) _F._modulus;
		}

		return res = y;
	}
}

template <class Vector1, class Vector2>
inline uint16 &DotProductDomain<Modular<uint16> >::dotSpecializedDD
	(uint16 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::const_iterator i = v1.begin ();
	typename Vector2::const_iterator j = v2.begin ();

	typename Vector1::const_iterator iterend = v1.begin () + v1.size() % _F._k;

	uint64 y = 0;

	for (; i != iterend; ++i, ++j)
		y += (uint64) *i * (uint64) *j;

	y %= (uint64) _F._modulus;

	for (; iterend != v1.end (); j += _F._k) {
		typename Vector1::const_iterator iter_i = iterend;
		typename Vector2::const_iterator iter_j;

		iterend += _F._k;

		for (iter_j = j; iter_i != iterend; ++iter_i, ++iter_j)
			y += (uint64) *iter_i * (uint64) *j;

		y %= (uint64) _F._modulus;
	}

	return res = y;
}

template <class Vector1, class Vector2>
inline uint16 &DotProductDomain<Modular<uint16> >::dotSpecializedDSP
	(uint16 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::first_type::const_iterator i_idx = v1.first.begin ();
	typename Vector1::second_type::const_iterator i_elt = v1.second.begin ();

	uint64 y = 0;

	if (v1.first.size () < _F._k) {
		for (; i_idx != v1.first.end (); ++i_idx, ++i_elt)
			y += (uint64) *i_elt * (uint64) v2[*i_idx];

		return res = y % (uint64) _F._modulus;
	} else {
		typename Vector1::first_type::const_iterator iterend = v1.first.begin () + v1.first.size() % _F._k;

		for (; i_idx != iterend; ++i_idx, ++i_elt)
			y += (uint64) *i_elt * (uint64) v2[*i_idx];

		y %= (uint64) _F._modulus;

		while (iterend != v1.first.end ()) {
			typename Vector1::first_type::const_iterator iter_i_idx = iterend;
			typename Vector1::second_type::const_iterator iter_i_elt = i_elt;

			iterend += _F._k;
			i_elt += _F._k;

			for (; iter_i_idx != iterend; ++iter_i_idx, ++iter_i_elt)
				y += (uint64) *iter_i_elt * (uint64) v2[*iter_i_idx];

			y %= (uint64) _F._modulus;
		}

		return res = y;
	}
}

template <class Vector1, class Vector2>
inline uint32 &DotProductDomain<Modular<uint32> >::dotSpecializedDD
	(uint32 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::const_iterator i;
	typename Vector2::const_iterator j;
  
	uint64 y = 0;
	uint64 t;

	for (i = v1.begin (), j = v2.begin (); i < v1.end (); ++i, ++j) {
		t = (uint64) *i * (uint64) *j;
		y += t;

		if (y < t)
			y += _F._two_64;
	}
  
	y %= (uint64) _F._modulus;

	return res = y;
}

template <class Vector1, class Vector2>
inline uint32 &DotProductDomain<Modular<uint32> >::dotSpecializedDSP
	(uint32 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::first_type::const_iterator i_idx;
	typename Vector1::second_type::const_iterator i_elt;
  
	uint64 y = 0;
	uint64 t;

	for (i_idx = v1.first.begin (), i_elt = v1.second.begin (); i_idx != v1.first.end (); ++i_idx, ++i_elt) {
		t = (uint64) *i_elt * (uint64) v2[*i_idx];
		y += t;

		if (y < t)
			y += _F._two_64;
	}
  
	y %= (uint64) _F._modulus;

	return res = y;
}

template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Modular<uint8> >::mulColDenseSpecialized
	(const VectorDomain<Modular<uint8> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::DenseVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::Column::const_iterator k;
	std::vector<uint32>::iterator l, l_end;

	if (_tmp.size () < w.size ())
		_tmp.resize (w.size ());

	std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

	l_end = _tmp.begin () + w.size ();

	do {
		j = v.begin ();
		j_end = j + __LINBOX_MIN (A->coldim (), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				*l += *k * *j;

		j_end += __LINBOX_MIN (A->coldim () - (j_end - v.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != v.end ());

	typename Vector1::iterator w_j;

	for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
		*w_j = *l;

	return w;
}

template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Modular<uint8> >::mulColDenseSpecialized
	(const VectorDomain<Modular<uint8> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseSequenceVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::Column::const_iterator k;
	std::vector<uint32>::iterator l, l_end;

	if (_tmp.size () < w.size ())
		_tmp.resize (w.size ());

	std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

	l_end = _tmp.begin () + w.size ();

	do {
		j = v.begin ();
		j_end = j + __LINBOX_MIN (A->coldim (), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				_tmp[k->first] += k->second * *j;

		j_end += __LINBOX_MIN (A->coldim () - (j_end - v.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != v.end ());

	typename Vector1::iterator w_j;

	for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
		*w_j = *l;

	return w;
}

template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Modular<uint8> >::mulColDenseSpecialized
	(const VectorDomain<Modular<uint8> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseAssociativeVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::Column::const_iterator k;
	std::vector<uint32>::iterator l, l_end;

	if (_tmp.size () < w.size ())
		_tmp.resize (w.size ());

	std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

	l_end = _tmp.begin () + w.size ();

	do {
		j = v.begin ();
		j_end = j + __LINBOX_MIN (A->coldim (), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				_tmp[k->first] += k->second * *j;

		j_end += __LINBOX_MIN (A->coldim () - (j_end - v.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != v.end ());

	typename Vector1::iterator w_j;

	for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
		*w_j = *l;

	return w;
}

template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Modular<uint8> >::mulColDenseSpecialized
	(const VectorDomain<Modular<uint8> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseParallelVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::Column::first_type::const_iterator k_idx;
	typename Matrix::Column::second_type::const_iterator k_elt;
	std::vector<uint32>::iterator l, l_end;

	if (_tmp.size () < w.size ())
		_tmp.resize (w.size ());

	std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

	l_end = _tmp.begin () + w.size ();

	do {
		j = v.begin ();
		j_end = j + __LINBOX_MIN (uint64 (A.coldim ()), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k_idx = i->first.begin (), k_elt = i->second.begin (), l = _tmp.begin ();
			     k_idx != i->first.end ();
			     ++k_idx, ++k_elt, ++l)
				_tmp[*k_idx] += *k_elt * *j;

		j_end += __LINBOX_MIN (uint64 (A.coldim () - (j_end - v.begin ())), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != v.end ());

	typename Vector1::iterator w_j;

	for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
		*w_j = *l;

	return w;
}

template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Modular<uint16> >::mulColDenseSpecialized
	(const VectorDomain<Modular<uint16> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::DenseVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j = v.begin (), j_end;
	typename Matrix::Column::const_iterator k;
	// Dan Roche, 7-1-04
	// std::vector<uint32>::iterator l, l_end;
	std::vector<uint64>::iterator l, l_end;

	if (_tmp.size () < w.size ())
		_tmp.resize (w.size ());

	std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

	l_end = _tmp.begin () + w.size ();

	do {
		j = v.begin ();
		j_end = j + __LINBOX_MIN (A->coldim (), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				*l += *k * *j;

		j_end += __LINBOX_MIN (A->coldim () - (j_end - v.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != v.end ());

	typename Vector1::iterator w_j;

	for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
		*w_j = *l;

	return w;
}

template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Modular<uint16> >::mulColDenseSpecialized
	(const VectorDomain<Modular<uint16> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseSequenceVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::Column::const_iterator k;
        // Dan Roche, 7-1-04
        // std::vector<uint32>::iterator l, l_end;
	std::vector<uint64>::iterator l, l_end;

	if (_tmp.size () < w.size ())
		_tmp.resize (w.size ());

	std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

	l_end = _tmp.begin () + w.size ();

	do {
		j = v.begin ();
		j_end = j + __LINBOX_MIN (A->coldim (), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				_tmp[k->first] += k->second * *j;

		j_end += __LINBOX_MIN (A->coldim () - (j_end - v.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != v.end ());

	typename Vector1::iterator w_j;

	for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
		*w_j = *l;

	return w;
}

template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Modular<uint16> >::mulColDenseSpecialized
	(const VectorDomain<Modular<uint16> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseAssociativeVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::Column::const_iterator k;
        // Dan Roche, 7-1-04
        // std::vector<uint32>::iterator l, l_end;
	std::vector<uint64>::iterator l, l_end;

	if (_tmp.size () < w.size ())
		_tmp.resize (w.size ());

	std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

	l_end = _tmp.begin () + w.size ();

	do {
		j = v.begin ();
		j_end = j + __LINBOX_MIN (A->coldim (), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
				_tmp[k->first] += k->second * *j;

		j_end += __LINBOX_MIN (A->coldim () - (j_end - v.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != v.end ());

	typename Vector1::iterator w_j;

	for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
		*w_j = *l;

	return w;
}

template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Modular<uint16> >::mulColDenseSpecialized
	(const VectorDomain<Modular<uint16> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseParallelVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j, j_end;
	typename Matrix::Column::first_type::const_iterator k_idx;
	typename Matrix::Column::second_type::const_iterator k_elt;
        // Dan Roche, 7-1-04
        // std::vector<uint32>::iterator l, l_end;
	std::vector<uint64>::iterator l, l_end;

	if (_tmp.size () < w.size ())
		_tmp.resize (w.size ());

	std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

	l_end = _tmp.begin () + w.size ();

	do {
		j = v.begin ();
		//Dan Roche, 7-2-04
		//j_end = j + __LINBOX_MIN (A->coldim (), VD.field ()._k);
		j_end = j + __LINBOX_MIN (A.coldim (), VD.field ()._k);

		for (; j != j_end; ++j, ++i)
			for (k_idx = i->first.begin (), k_elt = i->second.begin (), l = _tmp.begin ();
			     k_idx != i->first.end ();
			     ++k_idx, ++k_elt, ++l)
				_tmp[*k_idx] += *k_elt * *j;

		//j_end += __LINBOX_MIN (A->coldim () - (j_end - v.begin ()), VD.field ()._k);
		j_end += __LINBOX_MIN (A.coldim () - (j_end - v.begin ()), VD.field ()._k);

		for (l =_tmp.begin (); l != l_end; ++l)
			*l %= VD.field ()._modulus;

	} while (j_end != v.end ());

	typename Vector1::iterator w_j;

	for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
		*w_j = *l;

	return w;
}

template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Modular<uint32> >::mulColDenseSpecialized
	(const VectorDomain<Modular<uint32> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::DenseVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l;

	uint64 t;

	if (_tmp.size () < w.size ())
		_tmp.resize (w.size ());

	std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

	for (j = v.begin (); j != v.end (); ++j, ++i) {
		for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
			t = ((uint64) *k) * ((uint64) *j);

			*l += t;

			if (*l < t)
				*l += VD.field ()._two_64;
		}
	}

	typename Vector1::iterator w_j;

	for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
		*w_j = *l % VD.field ()._modulus;

	return w;
}

template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Modular<uint32> >::mulColDenseSpecialized
	(const VectorDomain<Modular<uint32> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseSequenceVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l;

	uint64 t;

	if (_tmp.size () < w.size ())
		_tmp.resize (w.size ());

	std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

	for (j = v.begin (); j != v.end (); ++j, ++i) {
		for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
			t = ((uint64) k->second) * ((uint64) *j);

			_tmp[k->first] += t;

			if (_tmp[k->first] < t)
				_tmp[k->first] += VD.field ()._two_64;
		}
	}

	typename Vector1::iterator w_j;

	for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
		*w_j = *l % VD.field ()._modulus;

	return w;
}

template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Modular<uint32> >::mulColDenseSpecialized
	(const VectorDomain<Modular<uint32> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseAssociativeVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j;
	typename Matrix::Column::const_iterator k;
	std::vector<uint64>::iterator l;

	uint64 t;

	if (_tmp.size () < w.size ())
		_tmp.resize (w.size ());

	std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

	for (j = v.begin (); j != v.end (); ++j, ++i) {
		for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
			t = ((uint64) k->second) * ((uint64) *j);

			_tmp[k->first] += t;

			if (_tmp[k->first] < t)
				_tmp[k->first] += VD.field ()._two_64;
		}
	}

	typename Vector1::iterator w_j;

	for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
		*w_j = *l % VD.field ()._modulus;

	return w;
}

template <class Vector1, class Matrix, class Vector2>
Vector1 &MVProductDomain<Modular<uint32> >::mulColDenseSpecialized
	(const VectorDomain<Modular<uint32> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseParallelVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j;
	typename Matrix::Column::first_type::const_iterator k_idx;
	typename Matrix::Column::second_type::const_iterator k_elt;
	std::vector<uint64>::iterator l;

	uint64 t;

	if (_tmp.size () < w.size ())
		_tmp.resize (w.size ());

	std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

	for (j = v.begin (); j != v.end (); ++j, ++i) {
		for (k_idx = i->first.begin (), k_elt = i->second.begin (), l = _tmp.begin ();
		     k_idx != i->first.end ();
		     ++k_idx, ++k_elt, ++l)
		{
			t = ((uint64) *k_elt) * ((uint64) *j);

			_tmp[*k_idx] += t;

			if (_tmp[*k_idx] < t)
				_tmp[*k_idx] += VD.field ()._two_64;
		}
	}

	typename Vector1::iterator w_j;

	for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
		*w_j = *l % VD.field ()._modulus;

	return w;
}

}

#endif // __LINBOX_field_modular_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
