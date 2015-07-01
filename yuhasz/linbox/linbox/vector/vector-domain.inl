/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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
 * See COPYING for license information.
 */

#ifndef __FIELD_VECTOR_DOMAIN_C
#define __FIELD_VECTOR_DOMAIN_C

#include "linbox-config.h"

#include <iostream>
#include <cctype>

#include "linbox/vector/vector-domain.h"
#include "linbox/util/debug.h"

namespace LinBox
{
	template <class Field>
	template <class Vector, class Trait>
	std::ostream &VectorDomain<Field>::writeSpecialized (std::ostream &os, const Vector &x,
							     VectorCategories::DenseVectorTag<Trait>) const
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
	std::ostream &VectorDomain<Field>::writeSpecialized (std::ostream &os, const Vector &x,
							     VectorCategories::SparseSequenceVectorTag<Trait>) const
	{
		typename Vector::const_iterator i;
		size_t idx;

		os << '[';

		for (i = x.begin (), idx = 0; i != x.end ();) {
			while (idx++ < i->first)
				os << "0, ";

			_F.write (os, i->second);

			if (++i != x.end ())
				os << ", ";
		}

		os << ']';

		return os;
	}

	template <class Field>
	template <class Vector, class Trait>
	std::ostream &VectorDomain<Field>::writeSpecialized (std::ostream &os, const Vector &x,
							     VectorCategories::SparseAssociativeVectorTag<Trait>) const
	{
		typename Vector::const_iterator i;
		size_t idx;

		os << '[';

		for (i = x.begin (), idx = 0; i != x.end ();) {
			while (idx++ < i->first)
				os << "0, ";

			_F.write (os, i->second);

			if (++i != x.end ())
				os << ", ";
		}

		os << ']';

		return os;
	}

	template <class Field>
	template <class Vector, class Trait>
	std::ostream &VectorDomain<Field>::writeSpecialized (std::ostream &os, const Vector &x,
							     VectorCategories::SparseParallelVectorTag<Trait>) const
	{
		typename Vector::first_type::const_iterator i;
		typename Vector::second_type::const_iterator j;
		size_t idx;

		os << '[';

		for (i = x.first.begin (), j = x.second.begin (), idx = 0; i != x.first.end ();) {
			while (idx++ < *i)
				os << "0, ";

			_F.write (os, *j);

			if (++i != x.first.end ())
				os << ", ";

			++j;
		}

		os << ']';

		return os;
	}

	template <class Field>
	template <class Vector, class Trait>
	std::istream &VectorDomain<Field>::readSpecialized (std::istream &is, const Vector &x,
							    VectorCategories::DenseVectorTag<Trait>) const
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
	std::istream &VectorDomain<Field>::readSpecialized (std::istream &is, const Vector &x,
							    VectorCategories::SparseSequenceVectorTag<Trait>) const
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
				x.push_back (std::pair <size_t, typename Field::Element> (idx, tmp));
			is >> c;
			idx++;
		}

		return is;
	}

	template <class Field>
	template <class Vector, class Trait>
	std::istream &VectorDomain<Field>::readSpecialized (std::istream &is, const Vector &x,
							    VectorCategories::SparseAssociativeVectorTag<Trait>) const
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
	template <class Vector, class Trait>
	std::istream &VectorDomain<Field>::readSpecialized (std::istream &is, const Vector &x,
							    VectorCategories::SparseParallelVectorTag<Trait>) const
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

			if (!_F.isZero (tmp)) {
				x.first.push_back (idx);
				x.second.push_back (tmp);
			}

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
			if (!VectorDomainBase<Field>::_F.areEqual (*i, *j))
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
		size_t idx;

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
		size_t idx;

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
						       VectorCategories::SparseParallelVectorTag<Trait1> tag1,
						       VectorCategories::DenseVectorTag<Trait2> tag2) const
	{
		typename Vector1::first_type::const_iterator i_idx;
		typename Vector1::second_type::const_iterator i_elt;
		typename Vector2::const_iterator j;
		size_t idx;

		for (i_idx = v1.first.begin (), i_elt = v1.second.begin (), j = v2.begin (), idx = 0;
		     i_idx != v1.first.end () && j != v2.end ();
		     j++, idx++)
		{
			if (*i_idx == idx) {
				if (!_F.areEqual (*i_elt, *j))
					return false;
				++i_idx;
				++i_elt;
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
						       VectorCategories::SparseParallelVectorTag<Trait1> tag1,
						       VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const
	{
		typename Vector1::first_type::const_iterator i_idx;
		typename Vector1::second_type::const_iterator i_elt;
		typename Vector2::const_iterator j;

		for (i_idx = v1.first.begin (), i_elt = v1.second.begin (), j = v2.begin ();
		     i_idx != v1.first.end () || j != v2.end ();)
		{
			while (i_idx != v1.first.end () && (j == v2.end () || *i_idx < j->first)) {
				if (!_F.isZero (*i_elt))
					return false;
				i_idx++;
				i_elt++;
			}

			while (j != v2.end () && (i_idx == v1.first.end () || j->first < *i_idx)) {
				if (!_F.isZero (j->second))
					return false;
				j++;
			}

			if (i_idx != v1.first.end () && j != v2.end () && *i_idx == j->first) {
				if (!_F.areEqual (*i_elt, j->second))
					return false;

				i_idx++; i_elt++; j++;
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
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	bool VectorDomain<Field>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						       VectorCategories::SparseParallelVectorTag<Trait1> tag1,
						       VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const
	{
		typename Vector1::first_type::const_iterator i_idx = v1.first.begin ();
		typename Vector1::second_type::const_iterator i_elt = v1.second.begin ();
		typename Vector2::const_iterator j = v2.begin ();

		while (i_idx != v1.first.end () || j != v2.end ()) {
			while (i_idx != v1.first.end () && (j == v2.end () || *i_idx < j->first)) {
				if (!_F.isZero (*i_elt))
					return false;
				i_idx++;
				i_elt++;
			}

			while (j != v2.end () && (i_idx == v1.first.end () || j->first < *i_idx)) {
				if (!_F.isZero (j->second))
					return false;
				j++;
			}

			if (i_idx != v1.first.end () && j != v2.end () && *i_idx == j->first) {
				if (!_F.areEqual (*i_elt, j->second))
					return false;

				i_idx++; i_elt++; j++;
			}
		}

		return true;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	bool VectorDomain<Field>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						       VectorCategories::SparseParallelVectorTag<Trait1> tag1,
						       VectorCategories::SparseParallelVectorTag<Trait2> tag2) const
	{
		typename Vector1::first_type::const_iterator i_idx = v1.first.begin ();
		typename Vector1::second_type::const_iterator i_elt = v1.second.begin ();
		typename Vector2::first_type::const_iterator j_idx = v2.first.begin ();
		typename Vector2::second_type::const_iterator j_elt = v2.second.begin ();

		while (i_idx != v1.first.end () || j_idx != v2.first.end ()) {
			while (i_idx != v1.first.end () && (j_idx == v2.first.end () || *i_idx < *j_idx)) {
				if (!_F.isZero (*i_elt))
					return false;
				i_idx++;
				i_elt++;
			}

			while (j_idx != v2.first.end () && (i_idx == v1.first.end () || *j_idx < *i_idx)) {
				if (!_F.isZero (*j_elt))
					return false;
				j_idx++;
				j_elt++;
			}

			if (i_idx != v1.first.end () && j_idx != v2.first.end () && *i_idx == *j_idx) {
				if (!_F.areEqual (*i_elt, *j_elt))
					return false;

				i_idx++; i_elt++; j_idx++; j_elt++;
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
	template <class Vector, class Trait>
	bool VectorDomain<Field>::isZeroSpecialized (const Vector &v, VectorCategories::SparseParallelVectorTag<Trait> tag) const
	{
		typename Vector::second_type::const_iterator i;

		for (i = v.second.begin (); i != v.second.end (); i++)
			if (!_F.isZero (*i))
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
				res.push_back (std::pair <size_t, typename Field::Element> (idx, *i));

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
						       VectorCategories::SparseParallelVectorTag<Trait1> tag1,
						       VectorCategories::DenseVectorTag<Trait2> tag2) const
	{
		typename Vector2::const_iterator i;
		int idx;

		res.first.clear ();
		res.second.clear ();

		for (i = v.begin (), idx = 0; i != v.end (); i++, idx++) {
			if (!_F.isZero (*i)) {
				res.first.push_back (idx);
				res.second.push_back (*i);
			}
		}

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
		size_t idx;

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
						       VectorCategories::SparseParallelVectorTag<Trait1> tag1,
						       VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const
	{
		typename Vector2::const_iterator i;

		res.first.clear ();
		res.second.clear ();

		for (i = v.begin (); i != v.end (); i++) {
			res.first.push_back (i->first);
			res.second.push_back (i->second);
		}

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
		size_t idx;

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
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v,
						       VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
						       VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const
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
						       VectorCategories::SparseParallelVectorTag<Trait1> tag1,
						       VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const
	{
		typename Vector2::const_iterator i;

		res.first.clear ();
		res.second.clear ();

		for (i = v.begin (); i != v.end (); i++) {
			res.first.push_back (i->first);
			res.second.push_back (i->second);
		}

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v,
						       VectorCategories::DenseVectorTag<Trait1> tag1,
						       VectorCategories::SparseParallelVectorTag<Trait2> tag2) const
	{
		typename Vector1::iterator i = res.begin ();
		typename Vector2::first_type::const_iterator j_idx = v.first.begin ();
		typename Vector2::second_type::const_iterator j_elt = v.second.begin ();
		size_t idx = 0;

		while (j_idx != v.first.end ()) {
			while (idx < *j_idx) {
				_F.init (*i, 0);
				++i; ++idx;
			}

			*i = *j_elt;

			++i; ++j_idx; ++j_elt; ++idx;
		}

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v,
						       VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
						       VectorCategories::SparseParallelVectorTag<Trait2> tag2) const
	{
		typename Vector2::first_type::const_iterator i_idx = v.first.begin ();
		typename Vector2::second_type::const_iterator i_elt = v.second.begin ();

		res.clear ();

		for (; i_idx != v.first.end (); ++i_idx, ++i_elt)
			res.push_back (std::pair <size_t, typename Field::Element> (*i_idx, *i_elt));

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v,
						       VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
						       VectorCategories::SparseParallelVectorTag<Trait2> tag2) const
	{
		typename Vector2::first_type::const_iterator i_idx = v.first.begin ();
		typename Vector2::second_type::const_iterator i_elt = v.second.begin ();

		res.clear ();

		for (; i_idx != v.first.end (); i_idx++, i_elt++)
			res[*i_idx] = *i_elt;

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v,
						       VectorCategories::SparseParallelVectorTag<Trait1> tag1,
						       VectorCategories::SparseParallelVectorTag<Trait2> tag2) const
	{
		res.first.resize (v.first.size ());
		std::copy (v.first.begin (), v.first.end (), res.first.begin ());
		res.second.resize (v.second.size ());
		std::copy (v.second.begin (), v.second.end (), res.second.begin ());
		return res;
	}

	template <class Field>
	template <class Vector1, class Trait, class Vector2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v, size_t i, size_t len,
						       VectorCategories::DenseVectorTag<Trait> tag) const
	{
		if (i == 0 && len == 0) {
			copy (res, v);
		} else {
			Vector1 res_part;

			copy (res_part, v);

			std::copy (res_part.begin (), (len == 0) ? res_part.end () : res_part.begin () + len, res.begin () + i);
		}

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait, class Vector2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v, size_t i, size_t len,
						       VectorCategories::SparseSequenceVectorTag<Trait> tag) const
	{
		if (i == 0 && len == 0)
			return copy (res, v);

		typename Vector1::iterator r_begin, r_end, part_end, iter;
		Vector1 res_part;

		copy (res_part, v);

		if (len == 0)
			part_end = res_part.end ();
		else
			part_end = std::lower_bound (res_part.begin (), res_part.end (), len,
						    VectorWrapper::CompareSparseEntries<Element> ());

		for (iter = res_part.begin (); iter != part_end; iter++)
			iter->first += i;

		r_begin = std::lower_bound (res.begin (), res.end (), i, VectorWrapper::CompareSparseEntries<Element> ());
		r_end = (len == 0) ? res.end () : std::lower_bound (r_begin, res.end (), i + len,
								    VectorWrapper::CompareSparseEntries<Element> ());
		r_begin = res.erase (r_begin, r_end);
		res.insert (r_begin, res_part.begin (), part_end);

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait, class Vector2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v, size_t i, size_t len,
						       VectorCategories::SparseAssociativeVectorTag<Trait> tag) const
	{
		if (i == 0 && len == 0)
			return copy (res, v);

		typename Vector1::iterator r_begin, r_end, part_end, iter;
		Vector1 res_part;

		copy (res_part, v);

		part_end = (len == 0) ? res_part.end () : res_part.find (len);

		r_begin = res.find (i);
		r_end = (len == 0) ? res.end () : res.find (i + len);
		res.erase (r_begin, r_end);

		for (iter = res_part.begin (); iter != part_end; iter++)
			res[iter->first + i] = iter->second;

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait, class Vector2>
	Vector1 &VectorDomain<Field>::copySpecialized (Vector1 &res, const Vector2 &v, size_t i, size_t len,
						       VectorCategories::SparseParallelVectorTag<Trait> tag) const
	{
		if (i == 0 && len == 0)
			return copy (res, v);

		typename Vector1::first_type::iterator r_idx_begin, r_idx_end, part_idx_end, iter;
		typename Vector1::second_type::iterator r_elt_begin, r_elt_end, part_elt_end;
		Vector1 res_part;

		copy (res_part, v);

		if (len == 0) {
			part_idx_end = res_part.first.end ();
			part_elt_end = res_part.second.end ();
		} else {
			part_idx_end = std::lower_bound (res_part.first.begin (), res_part.first.end (), len);
			part_elt_end = res_part.second.begin () + (part_idx_end - res_part.first.begin ());
		}

		for (iter = res_part.first.begin (); iter != part_idx_end; iter++)
			*iter += i;

		r_idx_begin = std::lower_bound (res.first.begin (), res.first.end (), i);
		r_elt_begin = res.second.begin () + (r_idx_begin - res.first.begin ());
		r_idx_end = (len == 0) ? res.first.end () : std::lower_bound (r_idx_begin, res.first.end (), i + len);
		r_idx_end = res.second.begin () + (r_idx_end - res.first.begin ());

		r_idx_begin = res.first.erase (r_idx_begin, r_idx_end);
		r_elt_begin = res.second.erase (r_elt_begin, r_elt_end);
		res.insert (r_idx_begin, res_part.first.begin (), part_idx_end);
		res.insert (r_elt_begin, res_part.second.begin (), part_elt_end);

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::copySpecialized (Vector &res, const Vector &v, size_t i, size_t len,
						      VectorCategories::DenseVectorTag<Trait> tag) const
	{
		std::copy (v.begin (), (len == 0) ? v.end () : v.begin () + len, res.begin () + i);
		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::copySpecialized (Vector &res, const Vector &v, size_t i, size_t len,
						      VectorCategories::SparseSequenceVectorTag<Trait> tag) const
	{
		typename Vector::const_iterator v_end;
		typename Vector::iterator r_begin, r_end, iter;
		typename Vector::difference_type offset;

		if (len == 0)
			v_end = v.end ();
		else
			v_end = std::lower_bound (v.begin (), v.end (), len,
						  VectorWrapper::CompareSparseEntries<Element> ());

		r_begin = std::lower_bound (res.begin (), res.end (), i, VectorWrapper::CompareSparseEntries<Element> ());
		r_end = (len == 0) ? res.end () : std::lower_bound (r_begin, res.end (), i + len,
								    VectorWrapper::CompareSparseEntries<Element> ());
		r_begin = res.erase (r_begin, r_end);
		offset = r_begin - res.begin ();
		res.insert (r_begin, v.begin (), v_end);
		r_begin = res.begin () + offset;
		r_end = r_begin + (v_end - v.begin ());

		for (iter = r_begin; iter != r_end; iter++)
			iter->first += i;

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::copySpecialized (Vector &res, const Vector &v, size_t i, size_t len,
						      VectorCategories::SparseAssociativeVectorTag<Trait> tag) const
	{
		typename Vector::const_iterator v_end;
		typename Vector::iterator r_begin, r_end, iter;

		v_end = (len == 0) ? v.end () : v.find (len);

		r_begin = res.find (i);
		r_end = (len == 0) ? res.end () : res.find (i + len);
		res.erase (r_begin, r_end);

		for (iter = v.begin (); iter != v_end; iter++)
			res[iter->first + i] = iter->second;

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::copySpecialized (Vector &res, const Vector &v, size_t i, size_t len,
						      VectorCategories::SparseParallelVectorTag<Trait> tag) const
	{
		typename Vector::first_type::const_iterator v_idx_end;
		typename Vector::second_type::const_iterator v_elt_end;
		typename Vector::first_type::iterator r_idx_begin, r_idx_end, iter;
		typename Vector::second_type::iterator r_elt_begin, r_elt_end;
		typename Vector::first_type::difference_type offset;

		if (len == 0) {
			v_idx_end = v.first.end ();
			v_elt_end = v.second.end ();
		} else {
			v_idx_end = std::lower_bound (v.first.begin (), v.first.end (), len);
			v_elt_end = v.second.begin () + (v_idx_end - v.first.begin ());
		}

		r_idx_begin = std::lower_bound (res.first.begin (), res.first.end (), i);
		r_elt_begin = res.second.begin () + (r_idx_begin - res.first.begin ());
		r_idx_end = (len == 0) ? res.first.end () : std::lower_bound (r_idx_begin, res.first.end (), i + len);
		r_elt_end = res.second.begin () + (r_idx_end - res.first.begin ());

		r_idx_begin = res.first.erase (r_idx_begin, r_idx_end);
		r_elt_begin = res.second.erase (r_elt_begin, r_elt_end);

		offset = r_idx_begin - res.first.begin ();
		res.first.insert (r_idx_begin, v.first.begin (), v_idx_end);
		res.second.insert (r_elt_begin, v.second.begin (), v_elt_end);

		r_idx_begin = res.first.begin () + offset;
		r_idx_end = r_idx_begin + (v_idx_end - v.first.begin ());

		for (iter = r_idx_begin; iter != r_idx_end; iter++)
			*iter += i;

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
	Vector1 &VectorDomain<Field>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						      VectorCategories::DenseVectorTag<Trait1>,
						      VectorCategories::DenseVectorTag<Trait2>,
						      VectorCategories::DenseVectorTag<Trait3>) const
	{
		typename Vector2::const_iterator i;
		typename Vector3::const_iterator j;
		typename Vector1::iterator k;

		linbox_check (y.size () == x.size ());
		linbox_check (res.size () == x.size ());

		for (i = y.begin (), j = x.begin (), k = res.begin (); i != y.end (); i++, j++, k++)
			_F.add (*k, *i, *j);

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
	Vector1 &VectorDomain<Field>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						      VectorCategories::SparseSequenceVectorTag<Trait1>,
						      VectorCategories::SparseSequenceVectorTag<Trait2>,
						      VectorCategories::SparseSequenceVectorTag<Trait3>) const
	{
		typename Vector2::const_iterator i;
		typename Vector3::const_iterator j;
		Element tmp;

		res.clear ();

		for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first < j->first) {
				res.push_back (*i);
				i++;
			}

			if (i != y.end () && i->first == j->first) {
				_F.add (tmp, i->second, j->second);
				if (!_F.isZero (tmp))
					res.push_back (std::pair <size_t, Element> (j->first, tmp));
				i++;
			} else {
				res.push_back (*j);
			}
		}

		while (i != y.end ()) {
			res.push_back (*i);
			i++;
		}

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
	Vector1 &VectorDomain<Field>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						      VectorCategories::SparseAssociativeVectorTag<Trait1>,
						      VectorCategories::SparseAssociativeVectorTag<Trait2>,
						      VectorCategories::SparseAssociativeVectorTag<Trait3>) const
	{
		typename Vector2::const_iterator i;
		typename Vector3::const_iterator j;
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
	template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
	Vector1 &VectorDomain<Field>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						      VectorCategories::SparseParallelVectorTag<Trait1>,
						      VectorCategories::SparseParallelVectorTag<Trait2>,
						      VectorCategories::SparseParallelVectorTag<Trait3>) const
	{
		typename Vector2::first_type::const_iterator i_idx = y.first.begin ();
		typename Vector3::first_type::const_iterator j_idx = x.first.begin ();
		typename Vector2::second_type::const_iterator i_elt = y.second.begin ();
		typename Vector3::second_type::const_iterator j_elt = x.second.begin ();
		Element tmp;

		res.first.clear ();
		res.second.clear ();

		for (; j_idx != x.first.end (); ++j_idx, ++j_elt) {
			while (i_idx != y.first.end () && *i_idx < *j_idx) {
				res.first.push_back (*i_idx);
				res.second.push_back (*i_elt);
				++i_idx; ++i_elt;
			}

			if (i_idx != y.first.end () && *i_idx == *j_idx) {
				_F.add (tmp, *i_elt, *j_elt);
				if (!_F.isZero (tmp)) {
					res.first.push_back (*j_idx);
					res.second.push_back (tmp);
				}
				++i_idx; ++i_elt;
			} else {
				res.first.push_back (*j_idx);
				res.second.push_back (*j_elt);
			}
		}

		while (i_idx != y.first.end ()) {
			res.first.push_back (*i_idx);
			res.second.push_back (*i_elt);
			++i_idx; ++i_elt;
		}

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::addinSpecialized (Vector1 &y, const Vector2 &x,
							VectorCategories::DenseVectorTag<Trait1>,
							VectorCategories::DenseVectorTag<Trait2>) const
	{
		typename Vector1::iterator i;
		typename Vector2::const_iterator j;

		linbox_check (y.size () == x.size ());

		for (i = y.begin (), j = x.begin (); i != y.end (); i++, j++)
			_F.addin (*i, *j);

		return y;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::addinSpecialized (Vector1 &y, const Vector2 &x,
							VectorCategories::SparseSequenceVectorTag<Trait1>,
							VectorCategories::SparseSequenceVectorTag<Trait2>) const
	{
		Vector1 res;

		add (res, y, x);
		copy (y, res);
		return y;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::addinSpecialized (Vector1 &y, const Vector2 &x,
							VectorCategories::SparseAssociativeVectorTag<Trait1>,
							VectorCategories::SparseAssociativeVectorTag<Trait2>) const
	{
		typename Vector1::iterator i;
		typename Vector2::const_iterator j;

		for (i = y.begin (), j = x.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first < j->first) i++;

			if (i != y.end () && i->first == j->first)
				_F.addin (i->second, j->second);
			else
				y[j->first] = j->second;
		}

		return y;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::addinSpecialized (Vector1 &y, const Vector2 &x,
							VectorCategories::SparseParallelVectorTag<Trait1>,
							VectorCategories::SparseParallelVectorTag<Trait2>) const
	{
		Vector1 res;

		add (res, y, x);
		copy (y, res);
		return y;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
	Vector1 &VectorDomain<Field>::subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						      VectorCategories::DenseVectorTag<Trait1>,
						      VectorCategories::DenseVectorTag<Trait2>,
						      VectorCategories::DenseVectorTag<Trait3>) const
	{
		typename Vector2::const_iterator i;
		typename Vector3::const_iterator j;
		typename Vector1::iterator k;

		linbox_check (y.size () == x.size ());
		linbox_check (res.size () == x.size ());

		for (i = y.begin (), j = x.begin (), k = res.begin (); i != y.end (); i++, j++, k++)
			VectorDomainBase<Field>::_F.sub (*k, *i, *j);

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
	Vector1 &VectorDomain<Field>::subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						      VectorCategories::SparseSequenceVectorTag<Trait1>,
						      VectorCategories::SparseSequenceVectorTag<Trait2>,
						      VectorCategories::SparseSequenceVectorTag<Trait3>) const
	{
		typename Vector2::const_iterator i;
		typename Vector3::const_iterator j;
		Element tmp;

		res.clear ();

		for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first < j->first) {
				res.push_back (*i);
				i++;
			}

			if (i != y.end () && i->first == j->first) {
				_F.sub (tmp, i->second, j->second);
				if (!_F.isZero (tmp))
					res.push_back (std::pair <size_t, Element> (j->first, tmp));
				i++;
			} else {
				res.push_back (std::pair <size_t, Element> (j->first, _F.neg (tmp, j->second)));
			}
		}

		while (i != y.end ()) {
			res.push_back (*i);
			i++;
		}

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
	Vector1 &VectorDomain<Field>::subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						      VectorCategories::SparseAssociativeVectorTag<Trait1>,
						      VectorCategories::SparseAssociativeVectorTag<Trait2>,
						      VectorCategories::SparseAssociativeVectorTag<Trait3>) const
	{
		typename Vector2::const_iterator i;
		typename Vector3::const_iterator j;
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
	template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
	Vector1 &VectorDomain<Field>::subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						      VectorCategories::SparseParallelVectorTag<Trait1>,
						      VectorCategories::SparseParallelVectorTag<Trait2>,
						      VectorCategories::SparseParallelVectorTag<Trait3>) const
	{
		typename Vector2::first_type::const_iterator i_idx = y.first.begin ();
		typename Vector3::first_type::const_iterator j_idx = x.first.begin ();
		typename Vector2::second_type::const_iterator i_elt = y.second.begin ();
		typename Vector3::second_type::const_iterator j_elt = x.second.begin ();
		Element tmp;

		res.first.clear ();
		res.second.clear ();

		for (; j_idx != x.first.end (); ++j_idx, ++j_elt) {
			while (i_idx != y.first.end () && *i_idx < *j_idx) {
				res.first.push_back (*i_idx);
				res.second.push_back (*i_elt);
				++i_idx; ++i_elt;
			}

			if (i_idx != y.first.end () && *i_idx == *j_idx) {
				_F.sub (tmp, *i_elt, *j_elt);
				if (!_F.isZero (tmp)) {
					res.first.push_back (*j_idx);
					res.second.push_back (tmp);
				}
				++i_idx; ++i_elt;
			} else {
				res.first.push_back (*j_idx);
				res.second.push_back (_F.neg (tmp, *j_elt));
			}
		}

		while (i_idx != y.first.end ()) {
			res.first.push_back (*i_idx);
			res.second.push_back (*i_elt);
			++i_idx; ++i_elt;
		}

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::subinSpecialized (Vector1 &y, const Vector2 &x,
							VectorCategories::DenseVectorTag<Trait1>,
							VectorCategories::DenseVectorTag<Trait2>) const
	{
		typename Vector1::iterator i;
		typename Vector2::const_iterator j;

		linbox_check (y.size () == x.size ());

		for (i = y.begin (), j = x.begin (); i != y.end (); i++, j++)
			_F.subin (*i, *j);

		return y;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::subinSpecialized (Vector1 &y, const Vector2 &x,
							VectorCategories::SparseSequenceVectorTag<Trait1>,
							VectorCategories::SparseSequenceVectorTag<Trait2>) const
	{
		Vector1 res;

		sub (res, y, x);
		copy (y, res);
		return y;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::subinSpecialized (Vector1 &y, const Vector2 &x,
							VectorCategories::SparseAssociativeVectorTag<Trait1>,
							VectorCategories::SparseAssociativeVectorTag<Trait2>) const
	{
		typename Vector1::iterator i;
		typename Vector2::const_iterator j;
		Element tmp;

		for (i = y.begin (), j = x.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first < j->first) i++;

			if (i != y.end () && i->first == j->first)
				_F.subin (i->second, j->second);
			else
				y[j->first] = _F.neg (tmp, j->second);
		}

		return y;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::subinSpecialized (Vector1 &y, const Vector2 &x,
							VectorCategories::SparseParallelVectorTag<Trait1>,
							VectorCategories::SparseParallelVectorTag<Trait2>) const
	{
		Vector1 res;

		sub (res, y, x);
		copy (y, res);
		return y;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::negSpecialized (Vector1 &res, const Vector2 &x,
						      VectorCategories::DenseVectorTag<Trait1>,
						      VectorCategories::DenseVectorTag<Trait2>) const
	{
		typename Vector2::const_iterator j;
		typename Vector1::iterator k;

		linbox_check (res.size () == x.size ());

		for (j = x.begin (), k = res.begin (); j != x.end (); ++j, ++k)
			_F.neg (*k, *j);

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::negSpecialized (Vector1 &res, const Vector2 &x,
						      VectorCategories::SparseSequenceVectorTag<Trait1>,
						      VectorCategories::SparseSequenceVectorTag<Trait2>) const
	{
		typename Vector2::const_iterator j;
		Element tmp;

		res.clear ();

		for (j = x.begin (); j != x.end (); ++j)
			res.push_back (std::pair <size_t, Element> (j->first, _F.neg (tmp, j->second)));

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::negSpecialized (Vector1 &res, const Vector2 &x,
						      VectorCategories::SparseAssociativeVectorTag<Trait1>,
						      VectorCategories::SparseAssociativeVectorTag<Trait2>) const
	{
		typename Vector2::const_iterator j;
		Element tmp;

		res.clear ();

		for (j = x.begin (); j != x.end (); ++j)
			res[j->first] = _F.neg (tmp, j->second);

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::negSpecialized (Vector1 &res, const Vector2 &x,
						      VectorCategories::SparseParallelVectorTag<Trait1>,
						      VectorCategories::SparseParallelVectorTag<Trait2>) const
	{
		typename Vector2::second_type::const_iterator j;
		Element tmp;

		res.first.resize (x.first.size ());
		res.second.clear ();

		std::copy (x.first.begin (), x.first.end (), res.first.begin ());

		for (j = x.second.begin (); j != x.second.end (); ++j)
			res.second.push_back (_F.neg (tmp, *j));

		return res;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::neginSpecialized (Vector &y,
						       VectorCategories::DenseVectorTag<Trait> tag) const
	{
		typename Vector::iterator i;

		for (i = y.begin (); i != y.end (); ++i)
			_F.negin (*i);

		return y;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::neginSpecialized (Vector &y,
						       VectorCategories::SparseSequenceVectorTag<Trait> tag) const
	{
		typename Vector::iterator i;

		for (i = y.begin (); i != y.end (); ++i)
			_F.negin (i->second);

		return y;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::neginSpecialized (Vector &y,
						       VectorCategories::SparseAssociativeVectorTag<Trait> tag) const
	{
		typename Vector::iterator i;

		for (i = y.begin (); i != y.end (); ++i)
			_F.negin (i->second);

		return y;
	}

	template <class Field>
	template <class Vector, class Trait>
	Vector &VectorDomain<Field>::neginSpecialized (Vector &y,
						       VectorCategories::SparseParallelVectorTag<Trait> tag) const
	{
		typename Vector::second_type::iterator i;

		for (i = y.second.begin (); i != y.second.end (); ++i)
			_F.negin (*i);

		return y;
	}

	template <class Field>
	template <class Vector1, class Vector2, class Trait>
	Vector1 &VectorDomain<Field>::mulSpecialized
		(Vector1                                 &res,
		 const Vector2                           &x,
		 const typename Field::Element           &a,
		 VectorCategories::DenseVectorTag<Trait>  tag) const
	{
		typename Vector2::const_iterator i;
		typename Vector1::iterator j;

		linbox_check (res.size () == x.size ());
	
		for (i = x.begin (), j = res.begin (); i != x.end (); ++i, ++j)
			_F.mul (*j, *i, a);

		return res;
	}

	template <class Field>
	template <class Vector1, class Vector2, class Trait>
	Vector1 &VectorDomain<Field>::mulSpecialized
		(Vector1                                          &res,
		 const Vector2                                    &x,
		 const typename Field::Element                    &a,
		 VectorCategories::SparseSequenceVectorTag<Trait>  tag) const
	{
		typename Vector2::const_iterator i;
		Element tmp;

		res.clear ();

		if (_F.isZero (a))
			return res;

		for (i = x.begin (); i != x.end (); i++)
			res.push_back (std::pair <size_t, Element> (i->first, _F.mul (tmp, i->second, a)));

		return res;
	}

	template <class Field>
	template <class Vector1, class Vector2, class Trait>
	Vector1 &VectorDomain<Field>::mulSpecialized
		(Vector1                                             &res,
		 const Vector2                                       &x,
		 const typename Field::Element                       &a,
		 VectorCategories::SparseAssociativeVectorTag<Trait>  tag) const
	{
		typename Vector2::const_iterator i;
		Element tmp;

		res.clear ();

		if (_F.isZero (a))
			return res;

		for (i = x.begin (); i != x.end (); i++)
			res[i->first] = _F.mul (tmp, i->second, a);

		return res;
	}

	template <class Field>
	template <class Vector1, class Vector2, class Trait>
	Vector1 &VectorDomain<Field>::mulSpecialized
		(Vector1                                          &res,
		 const Vector2                                    &x,
		 const typename Field::Element                    &a,
		 VectorCategories::SparseParallelVectorTag<Trait>  tag) const
	{
		typename Vector2::first_type::const_iterator i_idx;
		typename Vector2::second_type::const_iterator i_elt;
		Element tmp;

		res.first.clear ();
		res.second.clear ();

		if (_F.isZero (a))
			return res;

		for (i_idx = x.first.begin (); i_idx != x.first.end (); ++i_idx)
			res.first.push_back (*i_idx);

		for (i_elt = x.second.begin (); i_elt != x.second.end (); ++i_elt)
			res.second.push_back (_F.mul (tmp, *i_elt, a));

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

		for (i = x.begin (); i != x.end (); i++)
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

		for (i = x.begin (); i != x.end (); i++)
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
	Vector &VectorDomain<Field>::mulinSpecialized
		(Vector                                           &x,
		 const typename Field::Element                    &a,
		 VectorCategories::SparseParallelVectorTag<Trait>  tag) const
	{
		typename Vector::second_type::iterator i;

		if (_F.isZero (a)) {
			x.first.clear ();
			x.second.clear ();
			return x;
		}

		for (i = x.second.begin (); i != x.second.end (); i++)
			_F.mulin (*i, a);

		return x;
	}

	template <class Field>
	template <class Vector1, class Vector2, class Vector3, class Trait>
	Vector1 &VectorDomain<Field>::axpySpecialized
		(Vector1                                 &res,
		 const Vector2                           &y,
		 const typename Field::Element           &a,
		 const Vector3                           &x,
		 VectorCategories::DenseVectorTag<Trait>  tag) const
	{
		typename Vector2::const_iterator i;
		typename Vector3::const_iterator j;
		typename Vector1::iterator k;

		linbox_check (y.size () == x.size ());
		linbox_check (res.size () == x.size ());

		for (i = y.begin (), j = x.begin (), k = res.begin (); i != y.end (); i++, j++, k++)
			_F.axpy (*k, a, *j, *i);

		return res;
	}

	template <class Field>
	template <class Vector1, class Vector2, class Vector3, class Trait>
	Vector1 &VectorDomain<Field>::axpySpecialized
		(Vector1                                          &res,
		 const Vector2                                    &y,
		 const typename Field::Element                    &a,
		 const Vector3                                    &x,
		 VectorCategories::SparseSequenceVectorTag<Trait>  tag) const
	{
		typename Vector2::const_iterator i;
		typename Vector3::const_iterator j;
		Element tmp;

		res.clear ();

		for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first < j->first) {
				res.push_back (*i);
				i++;
			}

			if (i != y.end () && i->first == j->first) {
				_F.axpy (tmp, a, j->second, i->second);
				i++;
			}
			else
				_F.mul (tmp, a, j->second);

			if (!_F.isZero (tmp))
				res.push_back (std::pair <size_t, Element> (j->first, tmp));
		}

		while (i != y.end ()) {
			res.push_back (*i);
			i++;
		}

		return res;
	}

	template <class Field>
	template <class Vector1, class Vector2, class Vector3, class Trait>
	Vector1 &VectorDomain<Field>::axpySpecialized
		(Vector1                                             &res,
		 const Vector2                                       &y,
		 const typename Field::Element                       &a,
		 const Vector3                                       &x,
		 VectorCategories::SparseAssociativeVectorTag<Trait>  tag) const
	{
		typename Vector2::const_iterator i;
		typename Vector3::const_iterator j;

		res.clear ();

		for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
			while (i != y.end () && i->first < j->first) {
				res[i->first] = i->second;
				i++;
			}

			if (i != y.end () && i->first == j->first) {
				_F.axpy (res[j->first], a, j->second, i->second);
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
	template <class Vector1, class Vector2, class Vector3, class Trait>
	Vector1 &VectorDomain<Field>::axpySpecialized
		(Vector1                                          &res, 	
		 const Vector2                                    &y,
		 const typename Field::Element                    &a,
		 const Vector3                                    &x,
		 VectorCategories::SparseParallelVectorTag<Trait>  tag) const
	{
		typename Vector2::first_type::const_iterator i_idx = y.first.begin ();
		typename Vector3::first_type::const_iterator j_idx = x.first.begin ();
		typename Vector2::second_type::const_iterator i_elt = y.second.begin ();
		typename Vector3::second_type::const_iterator j_elt = x.second.begin ();
		Element tmp;

		res.first.clear ();
		res.second.clear ();

		for (; j_idx != x.first.end (); ++j_idx, ++j_elt) {
			while (i_idx != y.first.end () && *i_idx < *j_idx) {
				res.first.push_back (*i_idx);
				res.second.push_back (*i_elt);
				++i_idx; ++i_elt;
			}

			if (i_idx != y.first.end () && *i_idx == *j_idx) {
				_F.axpy (tmp, a, *j_elt, *i_elt);
				if (!_F.isZero (tmp)) {
					res.first.push_back (*j_idx);
					res.second.push_back (tmp);
				}
				++i_idx; ++i_elt;
			} else {
				res.first.push_back (*j_idx);
				res.second.push_back (_F.mul (tmp, *j_elt, a));
			}
		}

		while (i_idx != y.first.end ()) {
			res.first.push_back (*i_idx);
			res.second.push_back (*i_elt);
			++i_idx; ++i_elt;
		}

		return res;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::axpyinSpecialized
		(Vector1                                 &y,
		 const typename Field::Element           &a,
		 const Vector2                           &x,
		 VectorCategories::DenseVectorTag<Trait1> tag1,
		 VectorCategories::DenseVectorTag<Trait2> tag2) const
	{
		typename Vector1::iterator i;
		typename Vector2::const_iterator j;

		linbox_check (y.size () == x.size ());

		for (i = y.begin (), j = x.begin (); i != y.end (); ++i, ++j)
			_F.axpyin (*i, a, *j);

		return y;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::axpyinSpecialized
		(Vector1                                 &y,
		 const typename Field::Element           &a,
		 const Vector2                           &x,
		 VectorCategories::DenseVectorTag<Trait1>,
		 VectorCategories::SparseSequenceVectorTag<Trait2>) const
	{
		typename Vector2::const_iterator j;

		for (j = x.begin (); j != x.end (); ++j)
			_F.axpyin (y[j->first], a, j->second);

		return y;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::axpyinSpecialized
		(Vector1                                 &y,
		 const typename Field::Element           &a,
		 const Vector2                           &x,
		 VectorCategories::DenseVectorTag<Trait1>,
		 VectorCategories::SparseAssociativeVectorTag<Trait2>) const
	{
		typename Vector2::const_iterator j;

		for (j = x.begin (); j != x.end (); ++j)
			_F.axpyin (y[j->first], a, j->second);

		return y;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::axpyinSpecialized
		(Vector1                                 &y,
		 const typename Field::Element           &a,
		 const Vector2                           &x,
		 VectorCategories::DenseVectorTag<Trait1>,
		 VectorCategories::SparseParallelVectorTag<Trait2>) const
	{
		typename Vector2::first_type::const_iterator j_idx = x.first.begin ();
		typename Vector2::second_type::const_iterator j_elt = x.second.begin ();

		for (; j_idx != x.first.end (); ++j_idx, ++j_elt)
			_F.axpyin (y[*j_idx], a, *j_elt);

		return y;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::axpyinSpecialized
		(Vector1                                          &y,
		 const typename Field::Element                    &a,
		 const Vector2                                    &x,
		 VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
		 VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const
	{
		Vector1 res;

		axpy (res, a, x, y);
		copy (y, res);
		return y;
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &VectorDomain<Field>::axpyinSpecialized
		(Vector1                                             &y,
		 const typename Field::Element                       &a,
		 const Vector2                                       &x,
		 VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
		 VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const
	{
		typename Vector1::iterator i;
		typename Vector2::const_iterator j;
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
	Vector1 &VectorDomain<Field>::axpyinSpecialized
		(Vector1                                          &y,
		 const typename Field::Element                    &a,
		 const Vector2                                    &x,
		 VectorCategories::SparseParallelVectorTag<Trait1> tag1,
		 VectorCategories::SparseParallelVectorTag<Trait2> tag2) const
	{
		Vector1 res;

		axpy (res, a, x, y);
		copy (y, res);
		return y;
	}

	template <class Field>
	template <class Vector1, class Vector2>
	inline typename Field::Element &DotProductDomain<Field>::dotSpecializedDD
		(Element                                  &res,
		 const Vector1                            &v1,
		 const Vector2                            &v2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;
		accu.reset();

		linbox_check (v1.size () == v2.size ());

		for (i = v1.begin (), j = v2.begin (); i != v1.end (); i++, j++)
			accu.accumulate (*i, *j);

		return accu.get (res);
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline typename Field::Element &VectorDomain<Field>::dotSpecialized
		(Element                                           &res,
		 const Vector1                                     &v1,
		 const Vector2                                     &v2,
		 VectorCategories::SparseSequenceVectorTag<Trait1>  tag1,
		 VectorCategories::DenseVectorTag<Trait2>           tag2) const
	{
		typename Vector1::const_iterator i;
		
		accu.reset();

		for (i = v1.begin (); i != v1.end (); i++)
			accu.accumulate (i->second, v2[i->first]);

		return accu.get (res);
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline typename Field::Element &VectorDomain<Field>::dotSpecialized
		(Element                                              &res,
		 const Vector1                                        &v1,
		 const Vector2                                        &v2,
		 VectorCategories::SparseAssociativeVectorTag<Trait1>  tag1,
		 VectorCategories::DenseVectorTag<Trait2>              tag2) const
	{
		typename Vector1::const_iterator i;
		accu.reset();

		for (i = v1.begin (); i != v1.end (); i++)
			accu.accumulate (i->second, v2[i->first]);

		return accu.get (res);
	}

	template <class Field>
	template <class Vector1, class Vector2>
	inline typename Field::Element &DotProductDomain<Field>::dotSpecializedDSP
		(Element                                           &res,
		 const Vector1                                     &v1,
		 const Vector2                                     &v2) const
	{
		typename Vector1::first_type::const_iterator i_idx;
		typename Vector1::second_type::const_iterator i_elt;
		accu.reset();

		for (i_idx = v1.first.begin (), i_elt = v1.second.begin (); i_idx != v1.first.end (); ++i_idx, ++i_elt)
			accu.accumulate (*i_elt, v2[*i_idx]);

		return accu.get (res);
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline typename Field::Element &VectorDomain<Field>::dotSpecialized
		(Element                                              &res,
		 const Vector1                                        &v1,
		 const Vector2                                        &v2,
		 VectorCategories::SparseSequenceVectorTag<Trait1>     tag1,
		 VectorCategories::SparseSequenceVectorTag<Trait2>     tag2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;
		accu.reset();

		for (i = v1.begin (), j = v2.begin (); i != v1.end () && j != v2.end (); i++) {
			while (j != v2.end () && j->first < i->first) j++;

			if (j != v2.end () && j->first == i->first)
				accu.accumulate (i->second, j->second);
		}

		return accu.get (res);
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline typename Field::Element &VectorDomain<Field>::dotSpecialized
		(Element                                              &res,
		 const Vector1                                        &v1,
		 const Vector2                                        &v2,
		 VectorCategories::SparseAssociativeVectorTag<Trait1>  tag1,
		 VectorCategories::SparseSequenceVectorTag<Trait2>     tag2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;
		accu.reset();

		for (i = v1.begin (), j = v2.begin (); i != v1.end () && j != v2.end (); i++) {
			while (j != v2.end () && j->first < i->first) j++;

			if (j != v2.end () && j->first == i->first)
				accu.accumulate (i->second, j->second);
		}

		return accu.get (res);
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline typename Field::Element &VectorDomain<Field>::dotSpecialized
		(Element                                              &res,
		 const Vector1                                        &v1,
		 const Vector2                                        &v2,
		 VectorCategories::SparseParallelVectorTag<Trait1>     tag1,
		 VectorCategories::SparseSequenceVectorTag<Trait2>     tag2) const
	{
		typename Vector1::first_type::const_iterator i_idx = v1.first.begin ();
		typename Vector1::second_type::const_iterator i_elt = v1.second.begin ();
		typename Vector2::const_iterator j = v2.begin ();
		accu.reset();

		for (; i_idx != v1.first.end () && j != v2.end (); ++i_idx, ++i_elt) {
			while (j != v2.end () && j->first < *i_idx) j++;

			if (j != v2.end () && j->first == *i_idx)
				accu.accumulate (*i_elt, j->second);
		}

		return accu.get (res);
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline typename Field::Element &VectorDomain<Field>::dotSpecialized
		(Element                                              &res,
		 const Vector1                                        &v1,
		 const Vector2                                        &v2,
		 VectorCategories::SparseAssociativeVectorTag<Trait1>  tag1,
		 VectorCategories::SparseAssociativeVectorTag<Trait2>  tag2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;
		accu.reset();

		for (i = v1.begin (), j = v2.begin (); i != v1.end () && j != v2.end (); i++) {
			while (j != v2.end () && j->first < i->first) j++;

			if (j != v2.end () && j->first == i->first)
				accu.accumulate (i->second, j->second);
		}

		return accu.get (res);
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline typename Field::Element &VectorDomain<Field>::dotSpecialized
		(Element                                              &res,
		 const Vector1                                        &v1,
		 const Vector2                                        &v2,
		 VectorCategories::SparseParallelVectorTag<Trait1>     tag1,
		 VectorCategories::SparseAssociativeVectorTag<Trait2>  tag2) const
	{
		typename Vector1::first_type::const_iterator i_idx = v1.first.begin ();
		typename Vector1::second_type::const_iterator i_elt = v1.second.begin ();
		typename Vector2::const_iterator j = v2.begin ();
		accu.reset();

		for (; i_idx != v1.first.end () && j != v2.end (); ++i_idx, ++i_elt) {
			while (j != v2.end () && j->first < *i_idx) j++;

			if (j != v2.end () && j->first == *i_idx)
				accu.accumulate (*i_elt, j->second);
		}

		return accu.get (res);
	}

	template <class Field>
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline typename Field::Element &VectorDomain<Field>::dotSpecialized
		(Element                                              &res,
		 const Vector1                                        &v1,
		 const Vector2                                        &v2,
		 VectorCategories::SparseParallelVectorTag<Trait1>     tag1,
		 VectorCategories::SparseParallelVectorTag<Trait2>     tag2) const
	{
		typename Vector1::first_type::const_iterator i_idx = v1.first.begin ();
		typename Vector1::second_type::const_iterator i_elt = v1.second.begin ();
		typename Vector2::first_type::const_iterator j_idx = v2.first.begin ();
		typename Vector2::second_type::const_iterator j_elt = v2.second.begin ();
		accu.reset();

		for (; i_idx != v1.first.end () && j_idx != v2.first.end (); ++i_idx, ++i_elt) {
			while (j_idx != v2.first.end () && *j_idx < *i_idx) {
				j_idx++; j_elt++;
			}

			if (j_idx != v2.first.end () && *j_idx == *i_idx)
				accu.accumulate (*i_elt, *j_elt);
		}

		return accu.get (res);
	}

} // namespace LinBox

#endif // __FIELD_VECTOR_DOMAIN_H
