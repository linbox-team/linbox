/* linbox/matrix/dense-matrix.h
 * Copyright (C) 2014 the LinBox
 *
 * Written by :
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

/*! @file matrix/dense-matrix.h
 * @ingroup densematrix
 * @brief NODOC
 */

#ifndef __LINBOX_matrix_dense_matrix_H
#define __LINBOX_matrix_dense_matrix_H


//! @bug those are not just traits:
#include "linbox/linbox-config.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/matrix/matrix-traits.h"

namespace LinBox { /*  forward declarations */

	template <class _Field, class _Storage=typename Vector<_Field>::Dense >
	class BlasMatrix ;

	template <class _Field>
	using DenseMatrix = BlasMatrix<_Field> ;

	template<class _Matrix>
	class BlasSubmatrix ;

	template <class _Field>
	using DenseSubmatrix = BlasSubmatrix<DenseMatrix<_Field> > ;

	template <class _Matrix>
	class TriangularBlasMatrix ;

    template <class _Matrix>
	class TransposedBlasMatrix ;

    // template <class _Field, class _Storage=typename Vector<_Field>::Dense >
	// class TriangularBlasMatrix ;
    
	// template <class _Field>
	// using TriangularMatrix = TriangularBlasMatrix<_Field> ;

	/*! Write a matrix to a stream.
	 * The \c C++ way using <code>operator<< </code>
	 * @param o output stream
	 * @param Mat matrix to write.
	 */
	template <class _Field, class _Storage>
	std::ostream& operator<< (std::ostream & os, const BlasMatrix<_Field,_Storage> & Mat)
	{
		return Mat.write(os);
	}

	template <class _Matrix>
	std::ostream& operator<< (std::ostream & os, const BlasSubmatrix<_Matrix> & Mat)
	{
		return Mat.write(os);
	}

}

#include "linbox/matrix/densematrix/blas-matrix.h"
// #include "linbox/matrix/densematrix/blas-matrix-multimod.h"
// #include "linbox/matrix/densematrix/m4ri-matrix.h"

namespace LinBox { /*  MatrixContainerTrait */


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
	class MatrixContainerTrait<BlasSubmatrix<_Matrix> > {
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};

	template <class _Matrix>
	class MatrixContainerTrait<const BlasSubmatrix<_Matrix> > {
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};

} // LinBox

namespace LinBox { /*  MatrixTraits */

	template <class _Field, class _Rep>
	struct MatrixTraits< BlasMatrix<_Field,_Rep> > {
		typedef BlasMatrix<_Field,_Rep> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};

	template <class _Field, class _Rep>
	struct MatrixTraits< const BlasMatrix<_Field,_Rep> > {
		typedef const BlasMatrix<_Field,_Rep> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};

	template <class _Matrix>
	struct MatrixTraits< BlasSubmatrix<_Matrix> > {
		typedef BlasSubmatrix<_Matrix> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};

	template <class _Matrix>
	struct MatrixTraits< const BlasSubmatrix<_Matrix> > {
		typedef const BlasSubmatrix<_Matrix> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};

}

namespace LinBox { /*  MatrixHomTrait */

	//non working partial specialisation
	// template <class Ring, class Field, class _Rep>
	// struct MatrixHomTrait<BlasMatrix<Ring, _Rep>, Field> {
		// typedef BlasMatrix<Field,_Rep> value_type;
	// };

	template <class Ring, class Field>
	struct MatrixHomTrait<BlasMatrix<Ring, typename Vector<Ring>::Dense >, Field> {
		typedef BlasMatrix<Field,typename Vector<Field>::Dense > value_type;
	};


}

namespace LinBox { /*  IndexedCategory  */

	//! @bug this is trait, not a Category
	template<class Field, class _Rep>
	struct IndexedCategory< BlasMatrix<Field,_Rep> > {
		typedef IndexedTags::HasIndexed Tag;
	};

} // LinBox

namespace LinBox { /*  ContainerTraits */

	// this could also be a member of BlasVector
	template<class _Field, class _Rep>
	struct ContainerTraits<BlasMatrix<_Field,_Rep> > {
		typedef ContainerCategories::Matrix ContainerCategory ;
	};

	template<class _Matrix>
	struct ContainerTraits<BlasSubmatrix<_Matrix> > {
		typedef ContainerCategories::Matrix ContainerCategory ;
	};

	template<class _Matrix>
	struct ContainerTraits<TransposedBlasMatrix<_Matrix> > {
		//typedef ContainerCategories::Matrix ContainerCategory ;
        typedef typename ContainerTraits<_Matrix>::ContainerCategory ContainerCategory;
	};

}

#endif // __LINBOX_matrix_dense_matrix_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
