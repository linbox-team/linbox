/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/vector/vector-domain.inl
 * Copyright (C) 2001-2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-07-24 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Added support for the new SparseParallel vector type; this involves quite a
 * few new specializations.
 *
 * ------------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>
 *
 * Added the modifications for categories and vector traits that were designed
 * at the Rootbeer meeting. Added parametrization of VectorTags by VectorTraits.
 *
 * ------------------------------------
 * 2002-06-04 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Updated function definitions according to the new policy set in
 * vector-domain.h  This means the functions each take a parameter tag (or, for
 * dot product, tag1 and tag2) that allows specialization by vector type.
 * ------------------------------------
 *
 * 
 * ========LICENCE========
 * This file is part of the library LinBox.
 * 
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_field_vector_domain_gf2_INL
#define __LINBOX_field_vector_domain_gf2_INL


#include <iostream>
#include <cctype>

namespace LinBox
{ /*  VectorDomain */
	template <class Vector>
	std::ostream &VectorDomain<GF2>::writeSpecialized (std::ostream &os, const Vector &x,
							   VectorCategories::DenseZeroOneVectorTag) const
	{


		// TO BE REMOVED
		os << "writeSpec DenseZO, of size " << x.size() << ' ';

		os << "[ ";

		for (typename Vector::const_iterator i = x.begin (); i != x.end (); ++i)
			os << *i << ' ';

		os << ']';

		os << "( ";

		for (typename Vector::const_word_iterator i = x.wordBegin (); i != x.wordEnd (); ++i)
			os << *i << ' ';

		os << ')';

		return os;
	}

	template <class Vector>
	std::ostream &VectorDomain<GF2>::writeSpecialized (std::ostream &os, const Vector &x,
							   VectorCategories::SparseZeroOneVectorTag) const
	{
		typename Vector::const_iterator i;
		size_t idx = 0;

		// TO BE REMOVED
		os << "writeSpec SparseZO, of size " << x.size() << ' ';
		os << "[ ";

		for (i = x.begin (); i != x.end (); ++i) {
			while (++idx <= *i)
				os << 0 << ' ';

			os << 1 << ' ';
		}
		os << ']';

		return os;
	}

	template <class Vector>
	std::istream &VectorDomain<GF2>::readSpecialized (std::istream &is, const Vector &x,
							  VectorCategories::DenseZeroOneVectorTag) const
	{
		typename Vector::iterator i;
		char c;

		do { is >> c ; } while (!std::isdigit (c));

		is.unget ();

		for (i = x.begin (); i != x.end (); ++i)
			is >> *i;

		return is;
	}

	template <class Vector>
	std::istream &VectorDomain<GF2>::readSpecialized (std::istream &is, const Vector &x,
							  VectorCategories::SparseZeroOneVectorTag) const
	{
		char c;
		size_t idx;

		do { is >> c ; } while (!std::isdigit (c));

		is.unget ();
		x.clear ();

		while (1) {
			is >> c;

			if (!std::isdigit (c) && c != ' ') break;
			is.unget ();
			is >> idx;
			x.push_back (idx);
		}

		return is;
	}

	template <class Vector1, class Vector2>
	bool VectorDomain<GF2>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						     VectorCategories::DenseZeroOneVectorTag,
						     VectorCategories::SparseZeroOneVectorTag) const
	{
		typename Vector1::const_iterator i = v1.begin ();
		typename Vector2::const_iterator j = v2.begin ();
		size_t idx = 0;

		for (; j != v2.end (); ++j, ++i, ++idx) {
			while (idx < *j) {
				if (*i) return false;
				++idx;
				++i;
			}

			if (!*i) return false;
		}

		for (; i != v1.end (); ++i)
			if (*i) return false;

		return true;
	}

	template <class Vector1, class Vector2>
	bool VectorDomain<GF2>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						     VectorCategories::DenseZeroOneVectorTag,
						     VectorCategories::DenseZeroOneVectorTag) const
	{
		typename Vector1::const_word_iterator i = v1.wordBegin ();
		typename Vector2::const_word_iterator j = v2.wordBegin ();
		for (; j != v2.wordEnd (); ++j, ++i)
			if (*i != *j) return false;
		return true;
	}

	template <class Vector1, class Vector2>
	bool VectorDomain<GF2>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						     VectorCategories::SparseZeroOneVectorTag,
						     VectorCategories::SparseZeroOneVectorTag) const
	{
		return v1 == v2;
	}

	template <class Vector>
	bool VectorDomain<GF2>::isZeroSpecialized (const Vector &v,
						   VectorCategories::DenseZeroOneVectorTag) const
	{
		typename Vector::const_word_iterator i;

		for (i = v.wordBegin (); i != v.wordEnd (); ++i)
			if (*i) return false;

		return true;
	}

	template <class Vector1, class Vector2>
	Vector1 &VectorDomain<GF2>::copySpecialized (Vector1 &res, const Vector2 &v,
						     VectorCategories::SparseZeroOneVectorTag,
						     VectorCategories::DenseZeroOneVectorTag) const
	{
		typename Vector2::const_iterator i;
		size_t idx = 0;

		res.clear ();

		for (i = v.begin (); i != v.end (); ++i, ++idx)
			if (*i) res.push_back (idx);

		return res;
	}

	template <class Vector1, class Vector2>
	Vector1 &VectorDomain<GF2>::copySpecialized (Vector1 &res, const Vector2 &v,
						     VectorCategories::DenseZeroOneVectorTag,
						     VectorCategories::SparseZeroOneVectorTag) const
	{
		size_t sparsesize = *(v.rbegin());
		if (sparsesize > res.size()) res.resize( *(v.rbegin()) );
		std::fill (res.wordBegin (), res.wordEnd (), 0);

		for (typename Vector2::const_iterator i = v.begin ();
		     i != v.end ();
		     ++i)
			res[*i] = true;
		return res;
	}

	template <class Vector1, class Vector2>
	bool &VectorDomain<GF2>::dotSpecialized (bool &res, const Vector1 &v1, const Vector2 &v2,
						 VectorCategories::SparseZeroOneVectorTag,
						 VectorCategories::SparseZeroOneVectorTag) const
	{
		typename Vector1::const_iterator i = v1.begin ();
		typename Vector2::const_iterator j = v2.begin ();
		res = false;

		while (i != v1.end () || j != v2.end ()) {
			while (i != v1.end () && (j == v2.end () || *i < *j)) { res = !res; ++i; }
			while (j != v2.end () && (i == v1.end () || *j < *i)) { res = !res; ++j; }
			if (i != v1.end () && j != v2.end () && *i == *j) { ++i; ++j; }
		}

		return res;
	}

	template <class Vector1, class Vector2, class Vector3>
	Vector1 &VectorDomain<GF2>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						    VectorCategories::DenseZeroOneVectorTag,
						    VectorCategories::DenseZeroOneVectorTag,
						    VectorCategories::DenseZeroOneVectorTag) const
	{
		linbox_check (res.size () == y.size ());
		linbox_check (res.size () == x.size ());

		typename Vector1::word_iterator i = res.wordBegin ();
		typename Vector2::const_word_iterator j = y.wordBegin ();
		typename Vector3::const_word_iterator k = x.wordBegin ();

		for (; i != res.wordEnd (); ++i)
			*i = *j++ ^ *k++;

		return res;
	}

	template <class Vector1, class Vector2, class Vector3>
	Vector1 &VectorDomain<GF2>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						    VectorCategories::SparseZeroOneVectorTag,
						    VectorCategories::SparseZeroOneVectorTag,
						    VectorCategories::SparseZeroOneVectorTag) const
	{
		typename Vector2::const_iterator i = y.begin ();
		typename Vector3::const_iterator j = x.begin ();

		res.clear ();

		while (i != y.end () || j != x.end ()) {
			while (i != y.end () && (j == x.end () || *i < *j)) { res.push_back (*i); ++i; }
			while (j != x.end () && (i == y.end () || *j < *i)) { res.push_back (*j); ++j; }
			if (i != y.end () && j != x.end () && *i == *j) { ++i; ++j; }
		}

		return res;
	}

	template <class Vector1, class Vector2>
	Vector1 &VectorDomain<GF2>::addinSpecialized (Vector1 &y, const Vector2 &x,
						      VectorCategories::DenseZeroOneVectorTag,
						      VectorCategories::DenseZeroOneVectorTag) const
	{
		linbox_check (y.size () == x.size ());

		typename Vector1::word_iterator i = y.wordBegin ();
		typename Vector2::const_word_iterator j = x.wordBegin ();

		for (; i != y.wordEnd (); ++i, ++j)
			*i ^= *j;

		return y;
	}

	template <class Vector1, class Vector2>
	Vector1 &VectorDomain<GF2>::addinSpecialized (Vector1 &y, const Vector2 &x,
						      VectorCategories::DenseZeroOneVectorTag,
						      VectorCategories::SparseZeroOneVectorTag) const
	{
		typename Vector2::const_iterator i;

		for (i = x.begin (); i != x.end (); ++i)
			y[*i] = !y[*i];

		return y;
	}

}

namespace LinBox
{ /*  Dot product */
	// template<>
	template <class Vector1, class Vector2>
	inline bool &DotProductDomain<GF2>::dotSpecializedDD
	(bool          &res,
	 const Vector1 &v1,
	 const Vector2 &v2) const
	{
		linbox_check (v1.size () == v2.size ());

		unsigned long t = 0;
		typename Vector1::const_word_iterator i = v1.wordBegin ();
		typename Vector2::const_word_iterator j = v2.wordBegin ();

		while (i != v1.wordEnd () - 1)
			t ^= *i++ & *j++;

		const size_t zeroing = __LINBOX_BITSOF_LONG - (v1.size() % __LINBOX_BITSOF_LONG);
		unsigned long lastdot = *i & *j;
		lastdot <<= zeroing;
		lastdot >>= zeroing;

		t ^= lastdot;
		return res = __LINBOX_PARITY(t);
	}

	// template<>
	template <class Vector1, class Vector2>
	inline bool &DotProductDomain<GF2>::dotSpecializedDSP
	(bool          &res,
	 const Vector1 &v1,
	 const Vector2 &v2) const
	{
		typename Vector2::const_iterator i;

		res = 0;

		for (i = v2.begin (); i != v2.end (); ++i)
			res ^= v1[*i];

		return res;
	}

	// template<>
	template <class Vector1, class Vector2>
	inline BitVector::reference DotProductDomain<GF2>::dotSpecializedDD
	(BitVector::reference res,
	 const Vector1 &v1,
	 const Vector2 &v2) const
	{
		bool tmp;
		return res = dotSpecializedDD(tmp, v1, v2);
	}

	// template<>
	template <class Vector1, class Vector2>
	inline BitVector::reference DotProductDomain<GF2>::dotSpecializedDSP
	(BitVector::reference res,
	 const Vector1 &v1,
	 const Vector2 &v2) const
	{
		typename Vector2::const_iterator i;

		res = 0;

		for (i = v2.begin (); i != v2.end (); ++i)
			res ^= v1[*i];

		return res;
	}

}


#endif // __LINBOX_field_vector_domain_INL


