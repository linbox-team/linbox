/* linbox/matrix/sparse-coot-matrix.h
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

/*! @file matrix/SparseMatrix/sparse-coot-matrix.h
 * @ingroup matrix
 * @ingroup sparse
 * A <code>SparseMatrix2<_Field, SparseMatrixFormat::COO_T ></code> 
 * is a vector of (i,j,value) triples. 
 *
 * This is a variant of COO format.
 */

#ifndef __LINBOX_sparse_matrix_sparse_coot_matrix_H
#define __LINBOX_sparse_matrix_sparse_coot_matrix_H

#include "linbox/vector/vector-traits.h"
#include "linbox/blackbox/triplesbb.h"

namespace LinBox
{


	/** Sparse matrix, Coordinate storage.
	 *
	 * \ingroup matrix
	 * \ingroup sparse
	 */
	template<class _Field>
	class SparseMatrix2<_Field, SparseMatrixFormat::COO_T> 
	: public TriplesBB<_Field> {
	public :
		typedef _Field                             Field ; //!< Field
		typedef typename _Field::Element         Element ; //!< Element
		typedef const Element               constElement ; //!< const Element
		typedef SparseMatrixFormat::COO_T         Storage ; //!< Matrix Storage Format
		typedef SparseMatrix2<_Field,Storage>     Self_t ; //!< Self type
		typedef TriplesBB<Field> Father_t ; 
		SparseMatrix2() :  Father_t() {}

		SparseMatrix2(Field & F, size_t m, size_t n) 
		:  Father_t(F, m, n) {}

		void finalize(){}// this should be in Father_t

	};//sparse coot matrix

} // namespace LinBox

#endif // __LINBOX_sparse_matrix_sparse_coot_matrix_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
