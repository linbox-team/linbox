/* -*- mode: c; style: linux -*- */

/* linbox/field/vector-domain.C
 * Copyright (C) 2001-2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
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
	template <class Field, class Vector>
	Vector &VectorDomainBaseType (Dense)::mul
		(Vector                        &res,
		 const Vector                  &x,
		 const typename Field::Element &a) const
	{
		typename Vector::const_iterator i;
		typename Vector::iterator j;

		linbox_check (res.size () == x.size ());

		for (i = x.begin (), j = res.begin (); i < x.end (); i++, j++)
			_F.mul (*j, *i, a);

		return res;
	}

	template <class Field, class Vector>
	Vector &VectorDomainBaseType (SparseSequence)::mul
		(Vector                        &res,
		 const Vector                  &x,
		 const typename Field::Element &a) const
	{
		typename Vector::const_iterator i;
		typename Vector::iterator j;
		Element tmp;

		res.clear ();

		for (i = x.begin (); i < x.end (); i++)
			res.push_back (pair <size_t, Element> ((*i).first, _F.mul (tmp, (*i).second, a)));

		return res;
	}

	template <class Field, class Vector>
	Vector &VectorDomainBaseType (Dense)::mulin
		(Vector                        &x,
		 const typename Field::Element &a) const
	{
		typename Vector::iterator i;

		for (i = x.begin (); i < x.end (); i++)
			_F.mulin (*i, a);

		return x;
	}

	template <class Field, class Vector>
	Vector &VectorDomainBaseType (SparseSequence)::mulin
		(Vector                        &x,
		 const typename Field::Element &a) const
	{
		typename Vector::iterator i;

		for (i = x.begin (); i < x.end (); i++)
			_F.mulin ((*i).second, a);

		return x;
	}

	template <class Field, class Vector>
	Vector &VectorDomainBaseType (Dense)::axpy
		(Vector                        &res,
		 const Vector                  &y,
		 const typename Field::Element &a,
		 const Vector                  &x) const
	{
		typename Vector::const_iterator i, j;
		typename Vector::iterator k;

		linbox_check (y.size () == x.size ());
		linbox_check (res.size () == x.size ());

		for (i = y.begin (), j = x.begin (), k = res.begin (); i < y.end (); i++, j++, k++)
			_F.axpy (*k, *i, a, *j);

		return res;
	}

	template <class Field, class Vector>
	Vector &VectorDomainBaseType (Dense)::axpyin
		(Vector                        &y,
		 const typename Field::Element &a,
		 const Vector                  &x) const
	{
		typename Vector::iterator i;
		typename Vector::const_iterator j;

		linbox_check (y.size () == x.size ());

		for (i = y.begin (), j = x.begin (); i < y.end (); i++, j++)
			_F.axpyin (*i, a, *j);

		return y;
	}


	template <class Field, class Vector>
	Vector &VectorDomainBaseType (SparseSequence)::axpy
		(Vector                        &res,
		 const Vector                  &y,
		 const typename Field::Element &a,
		 const Vector                  &x) const
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

	template <class Field, class Vector>
	Vector &VectorDomainBaseType (SparseSequence)::axpyin
		(Vector                        &y,
		 const typename Field::Element &a,
		 const Vector                  &x) const
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

	template <class Field, class Vector1, class Vector2>
	typename VectorDomainType (Dense, Dense)::Element &VectorDomainType (Dense, Dense)::dotprod
		(Element       &res,
		 const Vector1 &v1,
		 const Vector2 &v2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;
		FieldAXPY<Field> r (_F);

		linbox_check (v1.size () == v2.size ());

		_F.init (res, 0);

		for (i = v1.begin (), j = v2.begin (); i < v1.end (); i++, j++)
			r.accumulate (*i, *j);

		return res = r.get ();
	}

	template <class Field, class Vector1, class Vector2>
	typename VectorDomainType (SparseSequence, Dense)::Element &VectorDomainType (SparseSequence, Dense)::dotprod
		(Element       &res,
		 const Vector1 &v1,
		 const Vector2 &v2) const
	{
		typename Vector1::const_iterator i;
		FieldAXPY<Field> r (_F);

		_F.init (res, 0);

		for (i = v1.begin (); i != v1.end (); i++)
			r.accumulate ((*i).second, v2[(*i).first]);

		return res = r.get ();
	}

} // namespace LinBox

#endif // __FIELD_VECTOR_DOMAIN_H
