/* -*- mode: c; style: linux -*- */

/* linbox/field/vector-domain.C
 * Copyright (C) 2001-2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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
 * See COPYING for license information.
 */

#ifndef __FIELD_VECTOR_DOMAIN_C
#define __FIELD_VECTOR_DOMAIN_C

#include "linbox-config.h"

#include <iostream>
#include <cctype>

#include "linbox/field/vector-domain.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"

namespace LinBox
{
	template <class Field>
	template <class Vector, class Trait>
	ostream &VectorDomain<Field>::writeSpecialized (ostream &os, const Vector &x,
							VectorCategories::DenseVectorTag<Trait> tag) const
	{
		typename Vector::const_iterator i;

		os << '[';

		for (i = x.begin (); i != x.end ();) {
			_F.write (os, *i);

			if (++i != x.end ())
				os << ", ";
		}

		os << ']';

		return os;
	}

	template <class Field>
	template <class Vector, class Trait>
	ostream &VectorDomain<Field>::writeSpecialized (ostream &os, const Vector &x,
							VectorCategories::SparseSequenceVectorTag<Trait> tag) const
	{
		typename Vector::const_iterator i;
		int idx;

		os << '[';

		for (i = x.begin (), idx = 0; i != x.end ();) {
			while (idx++ < i->first)
				os << "0, ";

			_F.write (os, *i);

			if (++i != x.end ())
				os << ", ";
		}

		os << ']';

		return os;
	}

	template <class Field>
	template <class Vector, class Trait>
	ostream &VectorDomain<Field>::writeSpecialized (ostream &os, const Vector &x,
							VectorCategories::SparseAssociativeVectorTag<Trait> tag) const
	{
		typename Vector::const_iterator i;
		int idx;

		os << '[';

		for (i = x.begin (), idx = 0; i != x.end ();) {
			while (idx++ < i->first)
				os << "0, ";

			_F.write (os, *i);

			if (++i != x.end ())
				os << ", ";
		}

		os << ']';

		return os;
	}

	template <class Field>
	template <class Vector, class Trait>
	istream &VectorDomain<Field>::readSpecialized (istream &is, const Vector &x,
						       VectorCategories::DenseVectorTag<Trait> tag) const
	{
		typename Vector::iterator i;
		char c;

		do is >> c; while (is && isspace (c));

		if (isdigit (c))
			is.unget (c);

		c = ',';

		i = x.begin ();

		while (is && c == ',') {
			do is >> c; while (is && isspace (c));
			is.unget (c);
			_F.read (is, *i++);
			is >> c;
		}

		return is;
	}

	template <class Field>
	template <class Vector, class Trait>
	istream &VectorDomain<Field>::readSpecialized (istream &is, const Vector &x,
						       VectorCategories::SparseSequenceVectorTag<Trait> tag) const
	{
		typename Field::Element tmp;
		char c;
		int idx;

		do is >> c; while (is && isspace (c));

		if (isdigit (c))
			is.unget (c);

		c = ','; x.clear (); idx = 0;

		while (is && c == ',') {
			do is >> c; while (is && isspace (c));
			is.unget (c);
			_F.read (is, tmp);
			if (!_F.isZero (tmp))
				x.push_back (pair <size_t, typename Field::Element> (idx, tmp));
			is >> c;
			idx++;
		}

		return is;
	}

	template <class Field>
	template <class Vector, class Trait>
	istream &VectorDomain<Field>::readSpecialized (istream &is, const Vector &x,
						       VectorCategories::SparseAssociativeVectorTag<Trait> tag) const
	{
		typename Field::Element tmp;
		char c;
		int idx;

		do is >> c; while (is && isspace (c));

		if (isdigit (c))
			is.unget (c);

		c = ','; x.clear (); idx = 0;

		while (is && c == ',') {
			do is >> c; while (is && isspace (c));
			is.unget (c);
			_F.read (is, tmp);
			if (!_F.isZero (tmp))
				x[idx] = tmp;
			is >> c;
			idx++;
		}

		return is;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	bool VectorDomain<Field>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						       VectorCategories::DenseVectorTag<Trait1> tag1,
						       VectorCategories::DenseVectorTag<Trait2> tag2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;

		if (v1.size () != v2.size ()) return false;

		for (i = v1.begin (), j = v2.begin (); i != v1.end (); i++, j++)
			if (!_F.areEqual (*i, *j))
				return false;

		return true;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	bool VectorDomain<Field>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						       VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
						       VectorCategories::DenseVectorTag<Trait2> tag2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;
		int idx;

		for (i = v1.begin (), j = v2.begin (), idx = 0; i != v1.end () && j != v2.end (); j++, idx++) {
			if (i->first == idx) {
				if (!_F.areEqual (i->second, *j))
					return false;
				i++;
			}
			else if (!_F.isZero (*j))
				return false;
		}

		return true;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	bool VectorDomain<Field>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						       VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
						       VectorCategories::DenseVectorTag<Trait2> tag2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;
		int idx;

		for (i = v1.begin (), j = v2.begin (), idx = 0; i != v1.end () && j != v2.end (); j++, idx++) {
			if (i->first == idx) {
				if (!_F.areEqual (i->second, *j))
					return false;
				i++;
			}
			else if (!_F.isZero (*j))
				return false;
		}

		return true;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	bool VectorDomain<Field>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						       VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
						       VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;

		for (i = v1.begin (), j = v2.begin (); i != v1.end () || j != v2.end ();) {
			while (i != v1.end () && (j == v2.end () || i->first < j->first)) {
				if (!_F.isZero (i->second))
					return false;
				i++;
			}

			while (j != v2.end () && (i == v1.end () || j->first < i->first)) {
				if (!_F.isZero (j->second))
					return false;
				j++;
			}

			if (i != v1.end () && j != v2.end () && i->first == j->first) {
				if (!_F.areEqual (i->second, j->second))
					return false;

				i++; j++;
			}
		}

		return true;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	bool VectorDomain<Field>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						       VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
						       VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;

		for (i = v1.begin (), j = v2.begin (); i != v1.end () || j != v2.end ();) {
			while (i != v1.end () && (j == v2.end () || i->first < j->first)) {
				if (!_F.isZero (i->second))
					return false;
				i++;
			}

			while (j != v2.end () && (i == v1.end () || j->first < i->first)) {
				if (!_F.isZero (j->second))
					return false;
				j++;
			}

			if (i != v1.end () && j != v2.end () && i->first == j->first) {
				if (!_F.areEqual (i->second, j->second))
					return false;

				i++; j++;
			}
		}

		return true;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	bool VectorDomain<Field>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						       VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
						       VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;

		for (i = v1.begin (), j = v2.begin (); i != v1.end () || j != v2.end ();) {
			while (i != v1.end () && (j == v2.end () || i->first < j->first)) {
				if (!_F.isZero (i->second))
					return false;
				i++;
			}

			while (j != v2.end () && (i == v1.end () || j->first < i->first)) {
				if (!_F.isZero (j->second))
					return false;
				j++;
			}

			if (i != v1.end () && j != v2.end () && i->first == j->first) {
				if (!_F.areEqual (i->second, j->second))
					return false;

				i++; j++;
			}
		}

		return true;
	}

	template <class Field>
	template <class Vector, class Trait>
	bool VectorDomain<Field>::isZeroSpecialized (const Vector &v, VectorCategories::DenseVectorTag<Trait> tag) const
	{
		typename Vector::const_iterator i;

		for (i = v.begin (); i != v.end (); i++)
			if (!_F.isZero (*i))
				return false;

		return true;
	}

	template <class Field>
	template <class Vector, class Trait>
	bool VectorDomain<Field>::isZeroSpecialized (const Vector &v, VectorCategories::SparseSequenceVectorTag<Trait> tag) const
	{
		typename Vector::const_iterator i;

		for (i = v.begin (); i != v.end (); i++)
			if (!_F.isZero (i->second))
				return false;

		return true;
	}

	template <class Field>
	template <class Vector, class Trait>
	bool VectorDomain<Field>::isZeroSpecialized (const Vector &v, VectorCategories::SparseAssociativeVectorTag<Trait> tag) const
	{
		typename Vector::const_iterator i;

		for (i = v.begin (); i != v.end (); i++)
			if (!_F.isZero (i->second))
				return false;

		return true;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v,
						       VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
						       VectorCategories::DenseVectorTag<Trait2> tag2) const
	{
		typename Vector2::const_iterator i;
		int idx;

		res.clear ();

		for (i = v.begin (), idx = 0; i != v.end (); i++, idx++)
			if (!_F.isZero (*i))
				res.push_back (pair <size_t, typename Field::Element> (idx, *i));

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v,
						       VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
						       VectorCategories::DenseVectorTag<Trait2> tag2) const
	{
		typename Vector2::const_iterator i;
		int idx;

		res.clear ();

		for (i = v.begin (), idx = 0; i != v.end (); i++, idx++)
			if (!_F.isZero (*i))
				res[idx] = *i;

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v,
						       VectorCategories::DenseVectorTag<Trait1> tag1,
						       VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const
	{
		typename Vector1::iterator i;
		typename Vector2::const_iterator j;
		int idx;

		for (i = res.begin (), j = v.begin (), idx = 0; j != v.end (); i++, j++, idx++) {
			while (idx < j->first) {
				_F.init (*i, 0);
				i++; idx++;
			}

			*i = j->second;
		}

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v,
						       VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
						       VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const
	{
		typename Vector2::const_iterator i;

		res.clear ();

		for (i = v.begin (); i != v.end (); i++)
			res[i->first] = i->second;

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v,
						       VectorCategories::DenseVectorTag<Trait1> tag1,
						       VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const
	{
		typename Vector1::iterator i;
		typename Vector2::const_iterator j;
		int idx;

		for (i = res.begin (), j = v.begin (), idx = 0; j != v.end (); i++, j++, idx++) {
			while (idx < j->first) {
				_F.init (*i, 0);
				i++; idx++;
			}

			*i = j->second;
		}

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v,
						       VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
						       VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const
	{
		typename Vector2::const_iterator i;

		res.clear ();

		for (i = v.begin (); i != v.end (); i++)
			res.push_back (*i);

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::addSpecialized (Vector &res, const Vector &y, const Vector &x,
						     VectorCategories::DenseVectorTag<Trait> tag) const
	{
		typename Vector::const_iterator i, j;
		typename Vector::iterator k;

		linbox_check (y.size () == x.size ());
		linbox_check (res.size () == x.size ());

		for (i = y.begin (), j = x.begin (), k = res.begin (); i < y.end (); i++, j++, k++)
			_F.add (*k, *i, *j);

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::addSpecialized (Vector &res, const Vector &y, const Vector &x,
						     VectorCategories::SparseSequenceVectorTag<Trait> tag) const
	{
		typename Vector::const_iterator i, j;
		Element tmp;

		res.clear ();

		for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first < j->first) {
				res.push_back (pair <size_t, Element> (i->first, i->second));
				i++;
			}

			if (i != y.end () && i->first == j->first) {
				_F.add (tmp, i->second, j->second);
				if (!_F.isZero (tmp))
					res.push_back (pair <size_t, Element> (j->first, tmp));
				i++;
			} else {
				res.push_back (pair <size_t, Element> (j->first, j->second));
			}
		}

		while (i != y.end ()) {
			res.push_back (pair <size_t, Element> (i->first, i->second));
			i++;
		}

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::addSpecialized (Vector &res, const Vector &y, const Vector &x,
						     VectorCategories::SparseAssociativeVectorTag<Trait> tag) const
	{
		typename Vector::const_iterator i, j;
		Element tmp;

		res.clear ();

		for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first < j->first) {
				res[i->first] = i->second;
				i++;
			}

			if (i != y.end () && i->first == j->first) {
				res[j->first] = _F.add (tmp, i->second, j->second);
				i++;
			} else {
				res[j->first] = j->second;
			}
		}

		while (i != y.end ()) {
			res[i->first] = i->second;
			i++;
		}

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::addinSpecialized (Vector &y, const Vector &x,
						       VectorCategories::DenseVectorTag<Trait> tag) const
	{
		typename Vector::iterator i;
		typename Vector::const_iterator j;

		linbox_check (y.size () == x.size ());

		for (i = y.begin (), j = x.begin (); i < y.end (); i++, j++)
			_F.addin (*i, *j);

		return y;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::addinSpecialized (Vector &y, const Vector &x,
						       VectorCategories::SparseSequenceVectorTag<Trait> tag) const
	{
		int i;
		typename Vector::const_iterator j;

		for (i = 0, j = x.begin (); j != x.end (); j++) {
			while (i < y.size () && y[i].first < j->first) i++;

			if (i < y.size () && y[i].first == j->first)
				_F.addin (y[i].second, j->second);
			else if (!_F.isZero (j->second))
				y.insert (y.begin () + i, *j);
		}

		return y;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::addinSpecialized (Vector &y, const Vector &x,
						       VectorCategories::SparseAssociativeVectorTag<Trait> tag) const
	{
		typename Vector::iterator i;
		typename Vector::const_iterator j;
		Element tmp;

		for (i = y.begin (), j = x.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first != j->first) i++;

			if (i != y.end () && i->first == j->first)
				_F.addin (i->second, j->second);
			else
				y[j->first] = j->second;
		}

		return y;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::subSpecialized (Vector &res, const Vector &y, const Vector &x,
						     VectorCategories::DenseVectorTag<Trait> tag) const
	{
		typename Vector::const_iterator i, j;
		typename Vector::iterator k;

		linbox_check (y.size () == x.size ());
		linbox_check (res.size () == x.size ());

		for (i = y.begin (), j = x.begin (), k = res.begin (); i < y.end (); i++, j++, k++)
			_F.sub (*k, *i, *j);

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::subSpecialized (Vector &res, const Vector &y, const Vector &x,
						     VectorCategories::SparseSequenceVectorTag<Trait> tag) const
	{
		typename Vector::const_iterator i, j;
		Element tmp;

		res.clear ();

		for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first < j->first) {
				res.push_back (pair <size_t, Element> (i->first, i->second));
				i++;
			}

			if (i != y.end () && i->first == j->first) {
				_F.sub (tmp, i->second, j->second);
				if (!_F.isZero (tmp))
					res.push_back (pair <size_t, Element> (j->first, tmp));
				i++;
			} else {
				res.push_back (pair <size_t, Element> (j->first, _F.neg (tmp, j->second)));
			}
		}

		while (i != y.end ()) {
			res.push_back (pair <size_t, Element> (i->first, i->second));
			i++;
		}

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::subSpecialized (Vector &res, const Vector &y, const Vector &x,
						     VectorCategories::SparseAssociativeVectorTag<Trait> tag) const
	{
		typename Vector::const_iterator i, j;
		Element tmp;

		res.clear ();

		for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first < j->first) {
				res[i->first] = i->second;
				i++;
			}

			if (i != y.end () && i->first == j->first) {
				res[j->first] = _F.sub (tmp, i->second, j->second);
				i++;
			} else {
				res[j->first] = _F.neg (tmp, j->second);
			}
		}

		while (i != y.end ()) {
			res[i->first] = i->second;
			i++;
		}

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::subinSpecialized (Vector &y, const Vector &x,
						       VectorCategories::DenseVectorTag<Trait> tag) const
	{
		typename Vector::iterator i;
		typename Vector::const_iterator j;

		linbox_check (y.size () == x.size ());

		for (i = y.begin (), j = x.begin (); i < y.end (); i++, j++)
			_F.subin (*i, *j);

		return y;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::subinSpecialized (Vector &y, const Vector &x,
						       VectorCategories::SparseSequenceVectorTag<Trait> tag) const
	{
		int i;
		typename Vector::const_iterator j;
		Element tmp;

		for (i = 0, j = x.begin (); j != x.end (); j++) {
			while (i < y.size () && y[i].first < j->first) i++;

			if (i < y.size () && y[i].first == j->first)
				_F.subin (y[i].second, j->second);
			else if (!_F.isZero (j->second))
				y.insert (y.begin () + i, pair <size_t, Element> (j->first, _F.neg (tmp, j->second)));
		}

		return y;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::subinSpecialized (Vector &y, const Vector &x,
						       VectorCategories::SparseAssociativeVectorTag<Trait> tag) const
	{
		typename Vector::iterator i;
		typename Vector::const_iterator j;
		Element tmp;

		for (i = y.begin (), j = x.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first != j->first) i++;

			if (i != y.end () && i->first == j->first)
				_F.subin (i->second, j->second);
			else
				y[j->first] = _F.neg (tmp, j->second);
		}

		return y;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::mulSpecialized
		(Vector                                  &res,
		 const Vector                            &x,
		 const typename Field::Element           &a,
		 VectorCategories::DenseVectorTag<Trait>  tag) const
	{
		typename Vector::const_iterator i;
		typename Vector::iterator j;

		linbox_check (res.size () == x.size ());
	
		for (i = x.begin (), j = res.begin (); i != x.end (); ++i, ++j)
			_F.mul (*j, *i, a);

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::mulSpecialized
		(Vector                                           &res,
		 const Vector                                     &x,
		 const typename Field::Element                    &a,
		 VectorCategories::SparseSequenceVectorTag<Trait>  tag) const
	{
		typename Vector::const_iterator i;
		typename Vector::iterator j;
		Element tmp;

		res.clear ();

		if (_F.isZero (a))
			return res;

		for (i = x.begin (); i < x.end (); i++)
			res.push_back (pair <size_t, Element> (i->first, _F.mul (tmp, i->second, a)));

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::mulSpecialized
		(Vector                                              &res,
		 const Vector                                        &x,
		 const typename Field::Element                       &a,
		 VectorCategories::SparseAssociativeVectorTag<Trait>  tag) const
	{
		typename Vector::const_iterator i;
		typename Vector::iterator j;
		Element tmp;

		res.clear ();

		if (_F.isZero (a))
			return res;

		for (i = x.begin (); i != x.end (); i++)
			res[i->first] = _F.mul (tmp, i->second, a);

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::mulinSpecialized
		(Vector                                  &x,
		 const typename Field::Element           &a,
		 VectorCategories::DenseVectorTag<Trait>  tag) const
	{
		typename Vector::iterator i;

		for (i = x.begin (); i < x.end (); i++)
			_F.mulin (*i, a);

		return x;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::mulinSpecialized
		(Vector                                           &x,
		 const typename Field::Element                    &a,
		 VectorCategories::SparseSequenceVectorTag<Trait>  tag) const
	{
		typename Vector::iterator i;

		if (_F.isZero (a)) {
			x.clear ();
			return x;
		}

		for (i = x.begin (); i < x.end (); i++)
			_F.mulin (i->second, a);

		return x;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::mulinSpecialized
		(Vector                                              &x,
		 const typename Field::Element                       &a,
		 VectorCategories::SparseAssociativeVectorTag<Trait>  tag) const
	{
		typename Vector::iterator i;

		if (_F.isZero (a)) {
			x.clear ();
			return x;
		}

		for (i = x.begin (); i != x.end (); i++)
			_F.mulin (i->second, a);

		return x;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::axpySpecialized
		(Vector                                  &res,
		 const Vector                            &y,
		 const typename Field::Element           &a,
		 const Vector                            &x,
		 VectorCategories::DenseVectorTag<Trait>  tag) const
	{
		typename Vector::const_iterator i, j;
		typename Vector::iterator k;

		linbox_check (y.size () == x.size ());
		linbox_check (res.size () == x.size ());

		for (i = y.begin (), j = x.begin (), k = res.begin (); i < y.end (); i++, j++, k++)
			_F.axpy (*k, *i, a, *j);

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::axpySpecialized
		(Vector                                           &res,
		 const Vector                                     &y,
		 const typename Field::Element                    &a,
		 const Vector                                     &x,
		 VectorCategories::SparseSequenceVectorTag<Trait>  tag) const
	{
		typename Vector::const_iterator i, j;
		Element tmp;

		res.clear ();

		for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first < j->first) {
				res.push_back (pair <size_t, Element> (i->first, i->second));
				i++;
			}

			if (i != y.end () && i->first == j->first) {
				_F.axpy (tmp, i->second, a, j->second);
				i++;
			}
			else
				_F.mul (tmp, a, j->second);

			if (!_F.isZero (tmp))
				res.push_back (pair <size_t, Element> (j->first, tmp));
		}

		while (i != y.end ()) {
			res.push_back (pair <size_t, Element> (i->first, i->second));
			i++;
		}

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::axpySpecialized
		(Vector                                              &res,
		 const Vector                                        &y,
		 const typename Field::Element                       &a,
		 const Vector                                        &x,
		 VectorCategories::SparseAssociativeVectorTag<Trait>  tag) const
	{
		typename Vector::const_iterator i, j;

		res.clear ();

		for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first < j->first) {
				res[i->first] = i->second;
				i++;
			}

			if (i != y.end () && i->first == j->first) {
				_F.axpy (res[j->first], i->second, a, j->second);
				i++;
			}
			else
				_F.mul (res[j->first], a, j->second);
		}

		while (i != y.end ()) {
			res[i->first] = i->second;
			i++;
		}

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::axpyinSpecialized
		(Vector                                  &y,
		 const typename Field::Element           &a,
		 const Vector                            &x,
		 VectorCategories::DenseVectorTag<Trait>  tag) const
	{
		typename Vector::iterator i;
		typename Vector::const_iterator j;

		linbox_check (y.size () == x.size ());

		for (i = y.begin (), j = x.begin (); i < y.end (); i++, j++)
			_F.axpyin (*i, a, *j);

		return y;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::axpyinSpecialized
		(Vector                                           &y,
		 const typename Field::Element                    &a,
		 const Vector                                     &x,
		 VectorCategories::SparseSequenceVectorTag<Trait>  tag) const
	{
		int i;
		typename Vector::const_iterator j;
		Element tmp;

		if (_F.isZero (a)) {
			y.clear ();
			return y;
		}

		for (i = 0, j = x.begin (); j != x.end (); j++) {
			while (i < y.size () && y[i].first < j->first) i++;

			if (i < y.size () && y[i].first == j->first)
				_F.axpyin (y[i].second, a, j->second);
			else if (!_F.isZero (j->second))
				y.insert (y.begin () + i, pair <size_t, Element> (j->first, _F.mul (tmp, a, j->second)));
		}

		return y;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::axpyinSpecialized
		(Vector                                              &y,
		 const typename Field::Element                       &a,
		 const Vector                                        &x,
		 VectorCategories::SparseAssociativeVectorTag<Trait>  tag) const
	{
		typename Vector::iterator i;
		typename Vector::const_iterator j;
		Element tmp;

		if (_F.isZero (a)) {
			y.clear ();
			return y;
		}

		for (i = y.begin (), j = x.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first != j->first) i++;

			if (i != y.end () && i->first == j->first)
				_F.axpyin (i->second, a, j->second);
			else if (!_F.isZero (j->second))
				y[j->first] = _F.mul (tmp, a, j->second);
		}

		return y;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	typename Field::Element &VectorDomain<Field>::dotSpecialized
		(Element                                  &res,
		 const Vector1                            &v1,
		 const Vector2                            &v2,
		 VectorCategories::DenseVectorTag<Trait1>  tag1,
		 VectorCategories::DenseVectorTag<Trait2>  tag2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;
		FieldAXPY<Field> r (_F);

		linbox_check (v1.size () == v2.size ());

		_F.init (res, 0);

		for (i = v1.begin (), j = v2.begin (); i < v1.end (); i++, j++)
			r.accumulate (*i, *j);

		return r.get (res);
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	typename Field::Element &VectorDomain<Field>::dotSpecialized
		(Element                                           &res,
		 const Vector1                                     &v1,
		 const Vector2                                     &v2,
		 VectorCategories::SparseSequenceVectorTag<Trait1>  tag1,
		 VectorCategories::DenseVectorTag<Trait2>           tag2) const
	{
		typename Vector1::const_iterator i;
		FieldAXPY<Field> r (_F);

		_F.init (res, 0);

		for (i = v1.begin (); i != v1.end (); i++)
			r.accumulate (i->second, v2[i->first]);

		return r.get (res);
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	typename Field::Element &VectorDomain<Field>::dotSpecialized
		(Element                                              &res,
		 const Vector1                                        &v1,
		 const Vector2                                        &v2,
		 VectorCategories::SparseAssociativeVectorTag<Trait1>  tag1,
		 VectorCategories::DenseVectorTag<Trait2>              tag2) const
	{
		typename Vector1::const_iterator i;
		FieldAXPY<Field> r (_F);

		_F.init (res, 0);

		for (i = v1.begin (); i != v1.end (); i++)
			r.accumulate (i->second, v2[i->first]);

		return r.get (res);
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	typename Field::Element &VectorDomain<Field>::dotSpecialized
		(Element                                              &res,
		 const Vector1                                        &v1,
		 const Vector2                                        &v2,
		 VectorCategories::SparseSequenceVectorTag<Trait1>     tag1,
		 VectorCategories::SparseSequenceVectorTag<Trait2>     tag2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;
		FieldAXPY<Field> r (_F);

		_F.init (res, 0);

		for (i = v1.begin (), j = v2.begin (); i != v1.end () && j != v2.end (); i++) {
			while (j != v2.end () && j->first < i->first) j++;

			if (j != v2.end () && j->first == i->first)
				r.accumulate (i->second, j->second);
		}

		return r.get (res);
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	typename Field::Element &VectorDomain<Field>::dotSpecialized
		(Element                                              &res,
		 const Vector1                                        &v1,
		 const Vector2                                        &v2,
		 VectorCategories::SparseAssociativeVectorTag<Trait1>  tag1,
		 VectorCategories::SparseSequenceVectorTag<Trait2>     tag2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;
		FieldAXPY<Field> r (_F);

		_F.init (res, 0);

		for (i = v1.begin (), j = v2.begin (); i != v1.end () && j != v2.end (); i++) {
			while (j != v2.end () && j->first < i->first) j++;

			if (j != v2.end () && j->first == i->first)
				r.accumulate (i->second, j->second);
		}

		return r.get (res);
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	typename Field::Element &VectorDomain<Field>::dotSpecialized
		(Element                                              &res,
		 const Vector1                                        &v1,
		 const Vector2                                        &v2,
		 VectorCategories::SparseAssociativeVectorTag<Trait1>  tag1,
		 VectorCategories::SparseAssociativeVectorTag<Trait2>  tag2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;
		FieldAXPY<Field> r (_F);

		_F.init (res, 0);

		for (i = v1.begin (), j = v2.begin (); i != v1.end () && j != v2.end (); i++) {
			while (j != v2.end () && j->first < i->first) j++;

			if (j != v2.end () && j->first == i->first)
				r.accumulate (i->second, j->second);
		}

		return r.get (res);
	}

} // namespace LinBox

#endif // __FIELD_VECTOR_DOMAIN_H
