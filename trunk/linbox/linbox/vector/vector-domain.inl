/* -*- mode: c; style: linux -*- */

/* linbox/field/matrix-domain.C
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

#ifndef __FIELD_MATRIX_DOMAIN_C
#define __FIELD_MATRIX_DOMAIN_C

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream>

#include "linbox/field/matrix-domain.h"

#include "linbox/debug.h"

namespace LinBox
{
	template <class Field, class Vector1, class Vector2>
	typename Field::element &MatrixDomainType (Dense, Dense)::dotprod
		(typename Field::element &res,
		 const Vector1           &v1,
		 const Vector2           &v2) const
	{
		typename Vector1::const_iterator i;
		typename Vector2::const_iterator j;
		element tmp;

		linbox_check (v1.size () == v2.size ());

		res = _F.zero ();

		for (i = v1.begin (), j = v2.begin (); i < v1.end (); i++, j++) {
			_F.mul (tmp, *i, *j);
			_F.addin (res, tmp);
		}

		return res;
	}

	template <class Field, class Vector1, class Vector2>
	typename Field::element &MatrixDomainType (SparseSequence, Dense)::dotprod
		(typename Field::element &res,
		 const Vector1           &v1,
		 const Vector2           &v2) const
	{
		typename Vector1::const_iterator i;
		element tmp;

		res = _F.zero ();

		for (i = v1.begin (); i < v1.end (); i++) {
			_F.mul (tmp, (*i).second, v2[(*i).first]);
			_F.addin (res, tmp);
		}

		return res;
	}

	template <class Field, class Vector1, class Vector2>
	Vector1 &MatrixDomainSimpleType (Dense)::axpy
		(Vector1                       &res,
		 const Vector1                 &y,
		 const typename Field::element &a,
		 const Vector1                 &x) const
	{
		typename Vector1::const_iterator i, j;
		typename Vector1::iterator k;
		element tmp;

		linbox_check (y.size () == x.size ());

		res.resize (y.size ());

		for (i = y.begin (), j = x.begin (), k = res.begin (); i < y.end (); i++, j++, k++) {
			_F.mul (tmp, a, *j);
			_F.add (*k, tmp, *i);
		}

		return res;
	}

	template <class Field, class Vector1, class Vector2>
	Vector1 &MatrixDomainSimpleType (Dense)::axpyin
		(Vector1                       &y,
		 const typename Field::element &a,
		 const Vector1                 &x) const
	{
		typename Vector1::iterator i;
		typename Vector1::const_iterator j;
		element tmp;

		linbox_check (y.size () == x.size ());

		res.resize (y.size ());

		for (i = y.begin (), j = x.begin (); i < y.end (); i++, j++) {
			_F.mul (tmp, a, *j);
			_F.addin (*i, tmp);
		}

		return y;
	}


	template <class Field, class Vector1, class Vector2>
	Vector1 &MatrixDomainSimpleType (SparseSequence)::axpy
		(Vector1                       &res,
		 const Vector1                 &y,
		 const typename Field::element &a,
		 const Vector1                 &x) const
	{
		typename Vector1::const_iterator i, j;
		typename Vector1::iterator k;
		element tmp;

		linbox_check (y.size () == x.size ());

		res.clear ();

		for (j = x.begin (), i = y.begin (); j < x.end (); j++) {
			mul (tmp, a, (*j).second);
			while ((*i).first < (*j).first)
				res.push_back (pair <size_t, element> ((*i).first, (*i).second));

			if ((*i).first == (*j).first)
				addin ((*i).second, tmp);
			else
				res.push_back (pair <size_t, element> ((*j).first, tmp));
		}

		return res;
	}

	template <class Field, class Vector1, class Vector2>
	Vector1 &MatrixDomainSimpleType (SparseSequence)::axpyin
		(Vector1                       &y,
		 const typename Field::element &a,
		 const Vector1                 &x) const
	{
		typename Vector1::iterator i;
		typename Vector1::const_iterator j;
		element tmp;

		linbox_check (y.size () == x.size ());

		res.resize (y.size ());

		for (i = y.begin (), j = x.begin (); j < x.end (); j++) {
			mul (tmp, a, (*j).second);
			while ((*i).first < (*j).first) i++;

			if ((*i).first == (*j).first)
				addin ((*i).second, tmp);
			else
				y.insert (i, pair <size_t, element> ((*j).first, tmp));
		}

		return y;
	}

} // namespace LinBox

#endif // __FIELD_MATRIX_DOMAIN_H
