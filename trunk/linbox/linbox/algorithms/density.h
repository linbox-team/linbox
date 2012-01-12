/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2010 LinBox
 *
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



#ifndef __LINBOX_density_H
#define __LINBOX_density_H
#include "linbox/vector/vector-traits.h"

namespace LinBox
{

	/** \brief Estimate nonzero entries in a vector, used in parallel elimination */
	template<class Vector>
	inline long density(const Vector& v)
	{

		return density(v, VectorTraits<Vector>::VectorCategory());
	}

	template<class Vector, class VectorCategory>
	inline long density(const Vector&, VectorCategory);

	template<class Vector>
	inline long density(const Vector& v, VectorCategories::DenseVectorTag)
	{

		return v.size();
	}


	template<class Vector>
	inline long density(const Vector& v, VectorCategories::SparseSequenceVectorTag)
	{

		return v.size();
	}


	template<class Vector>
	inline long density(const Vector& v, VectorCategories::SparseAssociativeVectorTag)
	{

		return v.size();
	}


	template<class Vector>
	inline long density(const Vector& v, VectorCategories::SparseParallelVectorTag) {


		return v.first.size();
	}
}

#endif //__LINBOX_density_H

