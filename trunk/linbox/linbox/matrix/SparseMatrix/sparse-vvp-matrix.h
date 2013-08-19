/* linbox/matrix/sparse-vvp-matrix.h
 * Copyright (C) 2013 the LinBox
 *
 * by -bds, templated from sparse-coo-matrix.h
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

/*! @file matrix/SparseMatrix/sparse-vvp-matrix.h
 * @ingroup matrix
 * @ingroup sparse
 * A <code>SparseMatrix2<_Field, SparseMatrixFormat::VVP ></code> 
 * is a vector or rows, 
 * wherein each row is a vector of index, (non-zero) value pairs.
 *
 * This is a variant of CSR format.
 */

#ifndef __LINBOX_sparse_matrix_sparse_vvp_matrix_H
#define __LINBOX_sparse_matrix_sparse_vvp_matrix_H

#include "linbox/vector/vector-traits.h"

namespace LinBox
{


	/** Sparse matrix, vector of row, with row = vector of index,value pairs.
	 * Variant of CSR storage.
	 *
	 * \ingroup matrix
	 * \ingroup sparse
	 */
	template<class _Field>
	class SparseMatrix2<_Field, SparseMatrixFormat::VVP> 
	: public SparseMatrix<_Field, typename Vector<_Field>::SparseSeq> {
	public :
		typedef _Field                             Field ; //!< Field
		typedef typename _Field::Element         Element ; //!< Element
		typedef const Element               constElement ; //!< const Element
		typedef SparseMatrixFormat::VVP         Storage ; //!< Matrix Storage Format
		typedef SparseMatrix2<_Field,Storage>     Self_t ; //!< Self type
		typedef SparseMatrix<_Field,typename Vector<_Field>::SparseSeq> Father_t ; 
		SparseMatrix2() :  Father_t() {}

		SparseMatrix2(Field & F, size_t m, size_t n) 
		:  Father_t(F, m, n) {}

		void finalize(){}// this should be in Father_t

	};//sparse vvp matrix

} // namespace LinBox

#endif // __LINBOX_sparse_matrix_sparse_vvp_matrix_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
