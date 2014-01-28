/* linbox/matrix/matrix-traits.h
 * Copyright (C) 2013 the LinBox group
 *
 * Written by :
 * bb <bbboyer@ncsu.edu>
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

#ifndef __LINBOX_matrix_matrix_traits_H
#define __LINBOX_matrix_matrix_traits_H

#include "linbox/util/debug.h"
#include "linbox/linbox-config.h"
#include "linbox/matrix/matrix-category.h"

namespace LinBox {

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

	//! @brief NO DOC.
	template <class Matrix>
	struct MatrixTraits {
		typedef Matrix MatrixType;
		typedef typename MatrixCategories::BlackboxTag MatrixCategory;
	};


} // LinBox

#endif // __LINBOX_matrix_matrix_traits_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
