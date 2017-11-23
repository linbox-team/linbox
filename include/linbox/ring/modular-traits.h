/* linbox/algorithms/sigma-basis.h
 * Copyright (C) 2017 Pascal Giorgi
 *
 * Written by Pascal Giorgi pascal.giorgi@lirmm.fr
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


#ifndef __LINBOX_modular_traits_h
#define __LINBOX_modular_traits_h


#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include <givaro/gfq.h>
#include "linbox/field/field-traits.h"


namespace LinBox
{ /* Givaro::Modular Base */

	using Givaro::Caster;

	template <class Ring>
	struct ClassifyRing;

	template <class Element, class Compute>
	struct ClassifyRing<Givaro::Modular<Element,Compute> const>
	{
		typedef RingCategories::ModularTag categoryTag;
	};

	template <class Element, class Compute>
	struct ClassifyRing<Givaro::Modular<Element,Compute>>
	{
		typedef RingCategories::ModularTag categoryTag;
	};

        template<class Element>
        struct ClassifyRing<Givaro::ModularBalanced<Element> > {
                typedef RingCategories::ModularTag categoryTag;
        };

        template<typename XXX>
        struct ClassifyRing<Givaro::GFqDom<XXX>> {
                typedef RingCategories::ModularTag categoryTag;
        };





/*
	template <>
	inline std::ostream& Givaro::ModularBase<Integer>::write (std::ostream &os) const
	{
		return os << "GMP integers mod " << _modulus;
	}
*/

	// template <>
	// inline integer& Givaro::Modular<integer>::init (integer& x, const double& y) const
	// {
	// 	integer tmp = (integer)y % _modulus;
	// 	if (tmp<0) tmp += _modulus;
	// 	return x = tmp;
	// }





} // namespace LinBox



#endif
