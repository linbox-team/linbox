/* linbox/matrix/matrix-category.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_matrix_category_H
#define __LINBOX_matrix_category_H


namespace LinBox
{

	template<class _Field, class _Rep>
	class BlasMatrix ;

	template<class _Matrix>
	class BlasSubmatrix ;

	template<class _Field, class _Storage>
	class SparseMatrix ;

	/** \brief For specializing matrix arithmetic
	 *
	 * This class defines matrix categories that allow us to specialize the matrix
	 * arithmetic in \ref MatrixDomain for different matrix representations. For
	 * example, a sparse matrix may have an efficient iterator over row vectors but
	 * not over column vectors. Therefore, an algorithm that tries to iterate over
	 * column vectors will run very slowly. Hence a specialization that avoids using
	 * column vectors is used instead.
	 */
	struct MatrixCategories {
		struct BlackboxTag { };
		struct RowMatrixTag : public virtual BlackboxTag { };
		struct ColMatrixTag : public virtual BlackboxTag { };
		struct RowColMatrixTag : public RowMatrixTag, public ColMatrixTag { };
	};


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

	template <class Field, class Rep>
	class MatrixContainerTrait<BlasMatrix<Field,Rep> > {
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};

	template <class Field, class Rep>
	class MatrixContainerTrait<const BlasMatrix<Field,Rep> > {
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};

	template <class _Matrix>
	class MatrixContainerTrait<const BlasSubmatrix<_Matrix> > {
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};

	template <class Field, class Storage>
	class MatrixContainerTrait<SparseMatrix<Field,Storage> > {
	public:
		typedef MatrixContainerCategory::Blackbox Type;
		// typedef MatrixContainerCategory::Container Type;
	};


}

#endif //__LINBOX_matrix_category_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
