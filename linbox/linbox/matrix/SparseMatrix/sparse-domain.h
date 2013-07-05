/* linbox/matrix/sparse-domain.h
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

/*! @file matrix/sparse-domain.h
 * @ingroup matrix
 * @ingroup sparse
 * A <code>SparseMatrix2<_Field ></code> ....
 */


#ifndef __LINBOX_sparse_matrix_sparse_domain_H
#define __LINBOX_sparse_matrix_sparse_domain_H

#include "linbox-config.h"
#include "linbox/util/debug.h"

namespace LinBox {

	//! @todo Vector knows Field
	template<class Field, class Vector>
	Vector & prepare(const Field & F , Vector & y, const typename Field::Element & a) {
		if ( !F.isOne(a) ) {
			if ( F.isZero(a) ) {
				for (size_t i = 0 ; i < y.size() ; ++i)
					F.assign(y[i],F.zero);
			}
			else if (F.isMOne(a)) {
				for (size_t i = 0 ; i < y.size() ; ++i)
					F.negin(y[i]);
			}
			else {
				for (size_t i = 0 ; i < y.size() ; ++i)
					F.mulin(y[i],a);
			}

		}
		return y ;
	}

} // LinBox

#endif // __LINBOX_sparse_matrix_sparse_domain_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
