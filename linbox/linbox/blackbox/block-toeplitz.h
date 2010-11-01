/* linbox/blackbox/block-toeplitz.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi pgiorgi@uwaterlo.ca
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_block_toeplitz_H
#define __LINBOX_block_toeplitz_H

#include <vector>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/algorithms/blas-domain.h>
#include <linbox/util/debug.h>
#include <linbox/blackbox/block-hankel.h>

namespace LinBox 
{
	
	template<class _Field>
	class BlockToeplitz : public BlockHankel<_Field> {
	public:
		typedef _Field Field;
		typedef typename Field::Element Element;
		
		BlockToeplitz(){}

		BlockToeplitz (const Field &F, const std::vector<BlasMatrix<Element> > &P, BlockHankelTag::shape s=BloackHankelTag::plain) 
			: BlockHankel(F, P, s) {}


		
	};

} // end of namespace LinBox

#endif

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
