/* linbox/matrix/sparse-domain.h
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

/*! @file matrix/sparsematrix/sparse-domain.h
 * @ingroup sparsematrix
 */


#ifndef __LINBOX_matrix_sparsematrix_sparse_domain_H
#define __LINBOX_matrix_sparsematrix_sparse_domain_H

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

namespace LinBox {

	/// y <- ay.  @todo Vector knows Field
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

#endif // __LINBOX_matrix_sparsematrix_sparse_domain_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
