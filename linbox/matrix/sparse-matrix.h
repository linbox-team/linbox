/* linbox/matrix/sparse-matrix.h
 * Copyright (C) 2013 the LinBox
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

/*! @file matrix/sparse-matrix.h
 * @ingroup matrix
 * @ingroup sparsematrix
 * A <code>SparseMatrix<_Field, _Storage></code> ....
 */

#ifndef __LINBOX_matrix_sparse_matrix_H
#define __LINBOX_matrix_sparse_matrix_H

#ifndef index_t
#define index_t ptrdiff_t
#endif

#include "linbox/matrix/sparse-formats.h"
#include "linbox/matrix/matrix-traits.h"

namespace LinBox {


	// Forward definition
	template<class _Field, class _Storage = SparseMatrixFormat::SparseSeq >
	class SparseMatrix ;


	template <class _Field, class _Storage>
	std::ostream& operator<< (std::ostream & os, const SparseMatrix<_Field,_Storage> & Mat)
	{
		return Mat.write(os);
	}

	template <class _Field, class _Storage>
	std::istream &operator >> (std::istream &is, SparseMatrix<_Field, _Storage> &A)
	{
		return A.read (is);
	}


	template<class _Field, class _Storage>
	class SparseMatrixDomain ;

} // LinBox

#include "linbox/matrix/sparsematrix/read-write-sparse.h"

#include "linbox/matrix/sparsematrix/sparse-generic.h"

#include "linbox/matrix/sparsematrix/sparse-coo-matrix.h"
// #include "linbox/matrix/sparsematrix/sparse-coo-1-matrix.h"
#include "linbox/matrix/sparsematrix/sparse-csr-matrix.h"
// #include "linbox/matrix/sparsematrix/sparse-csr-1-matrix.h"
#include "linbox/matrix/sparsematrix/sparse-ell-matrix.h"
#include "linbox/matrix/sparsematrix/sparse-ellr-matrix.h"
// #include "linbox/matrix/sparsematrix/sparse-ellr-1-matrix.h"
// #include "linbox/matrix/sparsematrix/sparse-bcsr-matrix.h"
// #include "linbox/matrix/sparsematrix/sparse-dia-matrix.h"
// #include "linbox/matrix/sparsematrix/sparse-hyb-matrix.h"
#include "linbox/matrix/sparsematrix/sparse-map-map-matrix.h"

#include "linbox/matrix/sparsematrix/sparse-tpl-matrix.h"
// #ifdef __LINBOX_USES_OPENMP
#ifdef _OPENMP
#include "linbox/matrix/sparsematrix/sparse-tpl-matrix-omp.h"
#endif

namespace LinBox { /*  MatrixContainerTraits */

	template <class Field, class Storage>
	class MatrixContainerTrait<SparseMatrix<Field,Storage> > {
	public:
                    //typedef MatrixContainerCategory::Blackbox Type;
                typedef MatrixContainerCategory::Container Type;
	};

	// template <class _Field, class _Storage>
	// struct MatrixTraits< SparseMatrix<_Field, _Storage> > {
		// typedef SparseMatrix<_Field, _Storage>      MatrixType;
		// typedef typename MatrixCategories::RowMatrixTag MatrixCategory;
	// };


	// template <class _Field, class _Storage>
	// struct GetEntryCategory<SparseMatrix<_Field,_Storage> > {
		  // typedef SolutionTags::Local Tag;
	// };
}

namespace LinBox { /*  MatrixTraits */

	template <class Field>
	struct MatrixTraits< SparseMatrix<Field, SparseMatrixFormat::CSR> >
	{
		typedef SparseMatrix<Field, SparseMatrixFormat::CSR> MatrixType;
		typedef typename MatrixCategories::IndexedMatrixTag MatrixCategory;
	};

	template <class Field>
	struct MatrixTraits< const SparseMatrix<Field, SparseMatrixFormat::CSR> >
	{
		typedef SparseMatrix<Field, SparseMatrixFormat::CSR> MatrixType;
		typedef typename MatrixCategories::IndexedMatrixTag MatrixCategory;
	};

	template <class Field>
	struct MatrixTraits< SparseMatrix<Field, SparseMatrixFormat::SparseMap> >
	{
		typedef SparseMatrix<Field, SparseMatrixFormat::SparseMap> MatrixType;
		typedef typename MatrixCategories::RowMatrixTag MatrixCategory;
	};

	template <class Field>
	struct MatrixTraits< const SparseMatrix<Field, SparseMatrixFormat::SparseMap> >
	{
		typedef SparseMatrix<Field, SparseMatrixFormat::SparseMap> MatrixType;
		typedef typename MatrixCategories::RowMatrixTag MatrixCategory;
	};

	template <class Field>
	struct MatrixTraits< SparseMatrix<Field, SparseMatrixFormat::SparsePar> >
	{
		typedef SparseMatrix<Field, SparseMatrixFormat::SparsePar> MatrixType;
		typedef typename MatrixCategories::RowMatrixTag MatrixCategory;
	};

	template <class Field>
	struct MatrixTraits< const SparseMatrix<Field, SparseMatrixFormat::SparsePar> >
	{
		typedef SparseMatrix<Field, SparseMatrixFormat::SparsePar> MatrixType;
		typedef typename MatrixCategories::RowMatrixTag MatrixCategory;
	};

	template <class Field>
	struct MatrixTraits< SparseMatrix<Field, SparseMatrixFormat::SparseSeq> >
	{
		typedef SparseMatrix<Field, SparseMatrixFormat::SparseSeq> MatrixType;
		typedef typename MatrixCategories::RowMatrixTag MatrixCategory;
	};

	template <class Field>
	struct MatrixTraits< const SparseMatrix<Field, SparseMatrixFormat::SparseSeq> >
	{
		typedef SparseMatrix<Field, SparseMatrixFormat::SparseSeq> MatrixType;
		typedef typename MatrixCategories::RowMatrixTag MatrixCategory;
	};

} // LinBox

namespace LinBox { /*  MatrixHomTrait */
	template <class Ring, class Field>
	struct MatrixHomTrait<SparseMatrix<Ring, SparseMatrixFormat::SparseSeq>, Field> {
		typedef SparseMatrix<Field, SparseMatrixFormat::SparseSeq> value_type;
	};

	template <class Ring, class Field>
	struct MatrixHomTrait<SparseMatrix<Ring, SparseMatrixFormat::SparsePar>, Field> {
		typedef SparseMatrix<Field, SparseMatrixFormat::SparsePar> value_type;
	};

	template <class Ring, class Field>
	struct MatrixHomTrait<SparseMatrix<Ring, SparseMatrixFormat::SparseMap>, Field> {
		typedef SparseMatrix<Field, SparseMatrixFormat::SparseMap> value_type;
	};

} // LinBox

namespace LinBox { /*  IndexedCategory */

	template<class Field, class Row>
	struct IndexedCategory< SparseMatrix<Field,Row> > 	{
		typedef IndexedTags::HasIndexed Tag;
	};

#if 1
	template<class Field>
	struct IndexedCategory< SparseMatrix<Field,SparseMatrixFormat::CSR> > 	{
		typedef IndexedTags::HasNext Tag;
	};

	template<class Field>
	struct IndexedCategory< SparseMatrix<Field,SparseMatrixFormat::COO> > 	{
		typedef IndexedTags::HasNext Tag;
	};

	template<class Field>
	struct IndexedCategory< SparseMatrix<Field,SparseMatrixFormat::ELL> > 	{
		typedef IndexedTags::HasNext Tag;
	};

	template<class Field>
	struct IndexedCategory< SparseMatrix<Field,SparseMatrixFormat::ELL_R> > 	{
		typedef IndexedTags::HasNext Tag;
	};

#endif




} // LinBox

namespace LinBox { /*  ContainerTraits */

	// this could also be a member of BlasVector
	template<class _Field>
	struct ContainerTraits<SparseMatrix<_Field, SparseMatrixFormat::COO> > {
		typedef ContainerCategories::Matrix ContainerCategory ;
	};

	template<class _Field>
	struct ContainerTraits<SparseMatrix<_Field, SparseMatrixFormat::CSR> > {
		typedef ContainerCategories::Matrix ContainerCategory ;
	};

	template<class _Field>
	struct ContainerTraits<SparseMatrix<_Field, SparseMatrixFormat::ELL> > {
		typedef ContainerCategories::Matrix ContainerCategory ;
	};

	template<class _Field>
	struct ContainerTraits<SparseMatrix<_Field, SparseMatrixFormat::ELL_R> > {
		typedef ContainerCategories::Matrix ContainerCategory ;
	};

}

namespace LinBox { /* Junk */

#if 1 /* make correspond SparseMatrixFormat::XXXSeq and Vector<Field>::XXXSeq */

		template<class Field>
		struct SparseVectorTranslate<Field,SparseMatrixFormat::SparseSeq> {
			typedef typename Vector<Field>::SparseSeq other_t;
		};

		template<class Field>
		struct SparseVectorTranslate<Field,typename Vector<Field>::SparseSeq> {
			typedef SparseMatrixFormat::SparseSeq other_t;
		};


		template<class Field>
		struct SparseVectorTranslate<Field,SparseMatrixFormat::SparsePar> {
			typedef typename Vector<Field>::SparsePar other_t;
		};

		template<class Field>
		struct SparseVectorTranslate<Field,SparseMatrixFormat::SparseMap> {
			typedef typename Vector<Field>::SparseMap other_t;
		};

#endif

} // LinBox

#endif // __LINBOX_matrix_sparse_matrix_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
