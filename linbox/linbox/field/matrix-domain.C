/* -*- mode: c; style: linux -*- */

/* linbox/src/field/matrix-domain.C
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
	MatrixDomain<Field>::Element &MatrixDomain<Field>::
		dotprod<Vector1, Vector2, VectorCategories::DenseVectorTag, VectorCategories::DenseVectorTag>
		(MatrixDomain<Field>::Element &res, const Vector1 &v1, const Vector2 &v2) const
	{
		Vector1::const_iterator i;
		Vector2::const_iterator j;
		Element tmp;

		linbox_check (v1.size () == v2.size ());

		res = _F.zero ();

		for (i = v1.begin (), j = v2.begin (); i < v1.end (); i++, j++) {
			_F.mul (tmp, *i, *j);
			_F.addin (res, tmp);
		}

		return res;
	}

	template <class Field, class Vector1, class Vector2>
	MatrixDomain<Field>::Element &MatrixDomain<Field>::
		dotprod<Vector1, Vector2, VectorCategories::SparseSequenceVectorTag, VectorCategories::DenseVectorTag>
		(MatrixDomain<Field>::Element &res, const Vector1 &v1, const Vector2 &v2) const
	{
		Vector1::const_iterator i;
		Element tmp;

		res = _F.zero ();

		for (i = v1.begin (); i < v1.end (); i++) {
			_F.mul (tmp, (*i).second, v2[(*i).first]);
			_F.addin (res, tmp);
		}

		return res;
	}

	template <class Field, class Vector>
	Vector &MatrixDomain<Field>::axpy<Vector, VectorCategories::DenseVectorTag>
		(Vector &res, const Vector &y, const MatrixDomain<Field>::Element &a, const Vector &x) const
	{
		Vector::const_iterator i, j;
		Vector::iterator k;
		Element tmp;

		linbox_check (y.size () == x.size ());

		res.resize (y.size ());

		for (i = y.begin (), j = x.begin (), k = res.begin (); i < y.end (); i++, j++, k++) {
			_F.mul (tmp, a, *j);
			_F.add (*k, tmp, *i);
		}

		return res;
	}

	template <class Field, class Vector>
	Vector &MatrixDomain<Field>::axpyin<Vector, VectorCategories::DenseVectorTag>
		(Vector &y, const MatrixDomain<Field>::Element &a, const Vector &x) const
	{
		Vector::iterator i;
		Vector::const_iterator j;
		Element tmp;

		linbox_check (y.size () == x.size ());

		res.resize (y.size ());

		for (i = y.begin (), j = x.begin (); i < y.end (); i++, j++) {
			_F.mul (tmp, a, *j);
			_F.addin (*i, tmp);
		}

		return y;
	}

} // namespace LinBox

#endif // __FIELD_MATRIX_DOMAIN_H
