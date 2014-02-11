/* Copyright (C) 2014 the LinBox group
 *
 *
 * Written by :
 *          BB <bbboyer@ncsu.edu>
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

/** @file linbox/blackbox/blackbox.h
 * @brief File including all of LinBox blackboxes.
 */


#ifndef __LINBOX_blackbox_blackbox_H
#define __LINBOX_blackbox_blackbox_H

#include "linbox/matrix/matrix-traits.h"

#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/lambda-sparse.h"
// #include "linbox/blackbox/subrowmatrix.h"
#include "linbox/blackbox/polynomial.h"
#include "linbox/blackbox/scalar-matrix.h"

namespace LinBox {

	template<class Field, class Row>
	struct IndexedCategory< LambdaSparseMatrix<Field,Row> >         {
		typedef IndexedTags::HasIndexed Tag;
	};


#if 0
	template<class Matrix, class MatrixCategory>
	struct IndexedCategory< SubRowMatrix<Matrix,MatrixCategory> >   {
		typedef IndexedTags::HasIndexed Tag;
	};
#endif

} // LinBox

#endif // __LINBOX_blackbox_blackbox_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
