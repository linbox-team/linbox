/* linbox/blackbox/block-toeplitz.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi pgiorgi@uwaterlo.ca
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

#ifndef __LINBOX_block_toeplitz_H
#define __LINBOX_block_toeplitz_H

#include <vector>
#include "linbox/matrix/blas-matrix.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/util/debug.h"
#include "linbox/blackbox/block-hankel.h"

namespace LinBox
{

	template<class _Field>
	class BlockToeplitz : public BlockHankel<_Field> {
	public:
		typedef _Field Field;
		typedef typename Field::Element Element;

		BlockToeplitz(){}

		BlockToeplitz (const Field &F, const std::vector<BlasMatrix<Element> > &P, BlockHankelTag::shape s=BloackHankelTag::plain) :
			BlockHankel(F, P, s) {}



	};

} // end of namespace LinBox

#endif


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
