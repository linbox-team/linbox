/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
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
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __MATRIX_CATEGORY_H
#define __MATRIX_CATEGORY_H

#include <linbox/matrix/sparse.h>
#include <linbox/matrix/dense.h>

#include <linbox/blackbox/sparse.h>
#include <linbox/blackbox/dense.h>


namespace LinBox{

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
	class MatrixContainerTrait<DenseMatrixBase<typename Field::Element> > {
	public:
		typedef MatrixContainerCategory::Container Type;
	};
	
	template <class Field>
	class MatrixContainerTrait<SparseMatrixBase<typename Field::Element> > {
	public:
		typedef MatrixContainerCategory::Container Type;
	};

	
	template <class Field>
	class MatrixContainerTrait<DenseMatrix<Field> > {
	public:
		typedef MatrixContainerCategory::Container Type;
	};

	template <class Field>
	class MatrixContainerTrait<SparseMatrix<Field> > {
	public:
		typedef MatrixContainerCategory::Container Type;
	};

}

#endif
