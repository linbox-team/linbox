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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream>

#include "linbox/field/vector-domain.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"

namespace LinBox
{
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

		for (i = x.begin (), j = res.begin (); i < x.end (); i++, j++)
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

		for (i = x.begin (); i < x.end (); i++)
			res.push_back (pair <size_t, Element> ((*i).first, _F.mul (tmp, (*i).second, a)));

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

		for (i = x.begin (); i < x.end (); i++)
			_F.mulin ((*i).second, a);

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

		for (j = x.begin (), i = y.begin (); j < x.end (); j++) {
			_F.mul (tmp, a, (*j).second);

			while (i < y.end () && (*i).first < (*j).first) {
				res.push_back (pair <size_t, Element> ((*i).first, (*i).second));
				i++;
			}

			if (i < y.end () && (*i).first == (*j).first) {
				_F.addin (tmp, (*i).second);
				i++;
			}

			res.push_back (pair <size_t, Element> ((*j).first, tmp));
		}

		while (i < y.end ()) {
			res.push_back (pair <size_t, Element> ((*i).first, (*i).second));
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
		Element tmp;

		res.clear ();

		for (j = x.begin (), i = y.begin (); j != x.end (); j++) {
			_F.mul (tmp, a, (*j).second);

			while (i != y.end () && (*i).first < (*j).first) {
				res[(*i).first] = (*i).second;
				i++;
			}

			if (i != y.end () && (*i).first == (*j).first) {
				_F.addin (tmp, (*i).second);
				i++;
			}

			res[(*j).first] = tmp;
		}

		while (i != y.end ()) {
			res[(*i).first] = (*i).second;
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
		typename Vector::iterator i;
		typename Vector::const_iterator j;
		Element tmp;

		for (i = y.begin (), j = x.begin (); j < x.end (); j++) {
			_F.mul (tmp, a, (*j).second);
			while (i < y.end () && (*i).first < (*j).first) i++;

			if (i < y.end () && (*i).first == (*j).first)
				_F.addin ((*i).second, tmp);
			else
				y.insert (i, pair <size_t, Element> ((*j).first, tmp));
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

		for (i = y.begin (), j = x.begin (); j != x.end (); j++) {
			_F.mul (tmp, a, (*j).second);
			while (i != y.end () && (*i).first != (*j).first) i++;

			if (i != y.end () && (*i).first == (*j).first)
				_F.addin ((*i).second, tmp);
			else
				y[(*j).first] = tmp;
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
			r.accumulate ((*i).second, v2[(*i).first]);

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
			r.accumulate ((*i).second, v2[(*i).first]);

		return r.get (res);
	}

} // namespace LinBox

#endif // __FIELD_VECTOR_DOMAIN_H
