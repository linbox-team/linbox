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


	//! @brief NO DOC.
	template <class Matrix>
	struct MatrixTraits {
		typedef Matrix MatrixType;
		typedef typename MatrixCategories::BlackboxTag MatrixCategory;
	};


	//! NODOC
	template <class Matrix>
	class MatrixContainerTrait {
	public:
		typedef MatrixContainerCategory::Blackbox Type;
	};


	//! try to map a blackbox over a homorphic ring
	//! The most suitable type
	template <class Blackbox, class Field>
	struct MatrixHomTrait {
		// static_assert(false,"should not be instanciated");
		//typedef ... FBlackbox
		// donot know
		// typedef Blackbox value_type;
	};

	/// limited doc so far
	namespace IndexedTags
	{
		struct HasIndexed{};
		struct HasNext{};
		struct NoIndexed{};
	}

	/// Trait to show whether or not the BB class has a Indexed iterator
	template<class BB>
	struct IndexedCategory {
		typedef IndexedTags::NoIndexed Tag;
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
