/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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

#ifndef __FIELD_MODULAR_INL
#define __FIELD_MODULAR_INL

#include <iostream>

namespace LinBox {

template <class Vector1, class Vector2>
inline uint16 &DotProductDomain<Modular<uint16> >::dotSpecializedDD
	(uint16 &res, const Vector1 &v1, const Vector2 &v2) const
{
	typename Vector1::const_iterator i = v1.begin ();
	typename Vector1::const_iterator iter_i;
	typename Vector1::const_iterator iterend;

	typename Vector2::const_iterator j = v2.begin ();
	typename Vector2::const_iterator iter_j;

	uint64 y = 0;

	iterend = v1.begin () + v1.size() % _F._k;

	for (; i != iterend; ++i, ++j)
		y += (uint64) *i * (uint64) *j;

	y %= (uint64) _F._modulus;

	for (; iterend != v1.end (); j += _F._k) {
		iter_i = iterend;
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
	typename Vector1::first_type::const_iterator iter_i_idx;
	typename Vector1::first_type::const_iterator iterend;
	typename Vector1::second_type::const_iterator i_elt = v1.second.begin ();
	typename Vector1::second_type::const_iterator iter_i_elt;

	uint64 y = 0;

	iterend = v1.first.begin () + v1.first.size() % _F._k;

	for (; i_idx != iterend; ++i_idx, ++i_elt)
		y += (uint64) *i_elt * (uint64) v2[*i_idx];

	y %= (uint64) _F._modulus;

	while (iterend != v1.first.end ()) {
		iter_i_idx = iterend;
		iter_i_elt = i_elt;
		iterend += _F._k;
		i_elt += _F._k;

		for (; iter_i_idx != iterend; ++iter_i_idx, ++iter_i_elt)
			y += (uint64) *iter_i_elt * (uint64) v2[*iter_i_idx];

		y %= (uint64) _F._modulus;
	}

	return res = y;
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

}

#endif // __FIELD_MODULAR_INL

