/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/matrix/matrix-domain.h
 * Copyright (C) 2002 Zhendong Wan, Bradford Hovinen
 *
 * Written by Zhendong Wan <wan@mail.eecis.udel.edu>,
 *            Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * ------------------------------------------------------------
 * 2002-11-26  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Added detailed documentation, cleaned up the interface slightly, and added
 * support for matrix traits. Added read, write, neg, negin, axpy, and
 * matrix-vector and matrix-black box operations.
 * ------------------------------------------------------------
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
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_matrix_domain_H
#define __LINBOX_matrix_domain_H

#include <iostream>

#include "linbox/field/gf2.h"
#include "linbox/matrix/matrix-domain.h"

// Specialization of MatrixDomain for GF2
namespace LinBox
{
	/*! Specialization of MatrixDomain for GF2.
	 * @bug this is half done and makes MatrixDomain on GF2 hardly usable.
	 * @todo this is where m4ri will play.
	 */
	template <>
	class MatrixDomain<GF2> {
	public:
		MatrixDomain (const GF2 &F) :
			_VD (F)
		{}

		template <class Vector1, class Matrix, class Vector2>
		Vector1 &vectorMul (Vector1 &w, const Matrix &A, const Vector2 &v) const
		{
			return mulSpecialized (w, A, v, typename MatrixTraits<Matrix>::MatrixCategory ());
		}

		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
					 MatrixCategories::RowMatrixTag) const
		{
			return mulRowSpecialized (w, A, v, typename VectorTraits<Vector1>::VectorCategory ());
		}

		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
					 MatrixCategories::ColMatrixTag) const
		{
			return mulColSpecialized (w, A, v,
						  typename VectorTraits<Vector1>::VectorCategory (),
						  typename VectorTraits<Vector2>::VectorCategory ());
		}
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
					 MatrixCategories::RowColMatrixTag) const
		{
			return mulRowSpecialized (w, A, v, typename VectorTraits<Vector1>::VectorCategory ());
		}

		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
					    VectorCategories::DenseZeroOneVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
					    VectorCategories::SparseZeroOneVectorTag) const;

		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
					    VectorCategories::DenseZeroOneVectorTag,
					    VectorCategories::DenseZeroOneVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
					    VectorCategories::DenseZeroOneVectorTag,
					    VectorCategories::SparseZeroOneVectorTag) const;

		VectorDomain<GF2> _VD;
	};

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MatrixDomain<GF2>::mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						       VectorCategories::DenseZeroOneVectorTag) const
	{
		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());

		typename Matrix::ConstRowIterator i = A.rowBegin ();
		typename Vector1::iterator j = w.begin ();

		for (; j != w.end (); ++j, ++i)
			_VD.dot (*j, v, *i);

		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MatrixDomain<GF2>::mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						       VectorCategories::SparseZeroOneVectorTag) const
	{
		typename Matrix::ConstRowIterator i = A.rowBegin ();
		GF2::Element t;
		unsigned int idx = 0;

		w.clear ();

		for (; i != A.rowEnd (); ++i, ++idx) {
			_VD.dot (t, v, *i);

			if (t)
				w.push_back (t);
		}

		return w;
	}

	template <class Vector1, class Matrix, class Vector2 >
	Vector1 &MatrixDomain<GF2>::mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						       VectorCategories::DenseZeroOneVectorTag,
						       VectorCategories::DenseZeroOneVectorTag) const
	{
		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());

		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j = v.begin ();

		//!@bug what's happening here ?
		_VD.subin (w, w);

		for (; j != v.end (); ++j, ++i)
			_VD.axpyin (w, *j, *i);

		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MatrixDomain<GF2>::mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						       VectorCategories::DenseZeroOneVectorTag,
						       VectorCategories::SparseZeroOneVectorTag) const
	{
		linbox_check (A.rowdim () == w.size ());

		typename Vector2::const_iterator j = v.begin ();

		_VD.subin (w, w);

		for (; j != v.end (); ++j) {
			typename Matrix::ConstColIterator i = A.colBegin () + *j;
			_VD.axpyin (w, true, *i);
		}

		return w;
	}

}

#endif // __LINBOX_matrix_domain_H

