/* linbox/matrix/sparse-matrix.h
 * Copyright (C) 2013 the LinBox
 *
 * Written by :
 * BB <bbboyer@ncsu.edu>
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
 * @ingroup sparse
 * A <code>SparseMatrix2<_Field ></code> ....
 */

#ifndef __LINBOX_matrix_sparse_matrix_H
#define __LINBOX_matrix_sparse_matrix_H


#include "linbox/matrix/sparse-formats.h"
#include "linbox/matrix/matrix-category.h"

namespace LinBox {
#if 0
	namespace Exception {

		/** Exception class for invalid matrix input.
		 * @todo should be factorised for all matrix readers.
		 */

		class InvalidMatrixInput {};
	} // Exception
#endif


	// Forward definition
	template<class _Field, class _Storage = SparseMatrixFormat::CSR >
	class SparseMatrix2 ;


	template <class _Field, class _Storage>
	std::istream &operator >> (std::istream &is, SparseMatrix2<_Field, _Storage> &A)
	{
		return A.read (is);
	}

	// template <class _Field, class _Storage>
	// struct MatrixTraits< SparseMatrix2<_Field, _Storage> > {
		// typedef SparseMatrix2<_Field, _Storage>      MatrixType;
		// typedef typename MatrixCategories::RowMatrixTag MatrixCategory;
	// };


	// template <class _Field, class _Storage>
	// struct GetEntryCategory<SparseMatrix2<_Field,_Storage> > {
		  // typedef SolutionTags::Local Tag;
	// };

} // LinBox

#include "SparseMatrix/sparse-generic.h"

#include "SparseMatrix/sparse-coo-matrix.h"
// #include "SparseMatrix/sparse-coo-1-matrix.h"
#include "SparseMatrix/sparse-csr-matrix.h"
// #include "SparseMatrix/sparse-csr-1-matrix.h"
#include "SparseMatrix/sparse-ell-matrix.h"
#include "SparseMatrix/sparse-ellr-matrix.h"
// #include "SparseMatrix/sparse-ellr-1-matrix.h"
// #include "SparseMatrix/sparse-bcsr-matrix.h"
// #include "SparseMatrix/sparse-dia-matrix.h"
#include "SparseMatrix/sparse-hyb-matrix.h"

#include "SparseMatrix/sparse-tpl-matrix.h"
// #ifdef __LINBOX_USES_OPENMP
#ifdef _OPENMP
#include "SparseMatrix/sparse-tpl-matrix-omp.h"
#endif


#endif // __LINBOX_matrix_sparse_matrix_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
