/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/matrix/matrix-category.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __LINBOX_matrix_category_H
#define __LINBOX_matrix_category_H

#include "linbox/matrix/sparse.h"
// #include "linbox/matrix/blas-matrix.h"

#include "linbox/blackbox/sparse.h"


namespace LinBox
{

	template<class _Field>
	class BlasMatrix ;

	template<class _Field>
	class BlasSubmatrix ;


	struct MatrixContainerCategory {
		struct BlasContainer{};
		struct Container{};
		struct Blackbox{};
	};

	template <class Matrix>
	class MatrixContainerTrait {
	public:
		typedef MatrixContainerCategory::Blackbox Type;
	};

	template <class Field>
	class MatrixContainerTrait<BlasMatrix<Field> > {
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};

	template <class Field>
	class MatrixContainerTrait<const BlasMatrix<Field> > {
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};

#if 0
	template <class Field>
	class MatrixContainerTrait<BlasMatrix<Field> > {
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};
#endif

	template <class Field>
	class MatrixContainerTrait<const BlasSubmatrix<Field> > {
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};

	template <class Element>
	class MatrixContainerTrait<SparseMatrixBase<Element> > {
	public:
		typedef MatrixContainerCategory::Container Type;
	};

	template <class Field>
	class MatrixContainerTrait<SparseMatrix<Field> > {
	public:
		typedef MatrixContainerCategory::Container Type;
	};

}

#endif //__LINBOX_matrix_category_H

