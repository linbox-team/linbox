/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/field/field-traits.h
 * Copyright (C) June 2004 Dan Roche
 *  ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_field_traits_H
#define __LINBOX_field_traits_H

#include "linbox/integer.h"

// Namespace in which all LinBox library code resides
namespace LinBox {

	/*! \brief some basic information about each field or ring.
	  \ingroup field

	 * It will try
	 * to take the information from the field when possible, and use defaults
	 * otherwise.
	 * maxModulus returns the greatest modulus that is usable for a given field, -1
	 * for infinite. goodModulus returns takes an integer and returns true if that
	 * integer is a valid modulus for the given field.
	 * maxExponent and goodExponent do the same for the prime power.
	 */
	namespace RingCategories {
		//! generic ring.
		struct GenericTag{};
		//! If it is isomorphic to Z/mZ, for some m or its extensions.
		struct ModularTag : public virtual GenericTag{};
		//! If it is isomorphic to Z
		struct IntegerTag : public virtual GenericTag{};
		//! If it is isomorphic to Q
		struct RationalTag : public virtual GenericTag{};
	}

	/*! Default ring category.
	 */
	template <class Field>
	struct ClassifyRing {
		typedef	RingCategories::GenericTag categoryTag;
	};

	/*! FieldTrait.
	 */
	template <class Field>
	struct FieldTraits {
		typedef typename ClassifyRing<Field>::categoryTag categoryTag;

		static integer& maxModulus( integer& i )
		{
			return i = static_cast<integer>(Field::getMaxModulus());
		}

		static uint64_t & maxModulus( uint64_t& i )
		{
			return i = static_cast<uint64_t>(Field::getMaxModulus());
		}

		static uint32_t & maxModulus( uint32_t& i )
		{
			return i = static_cast<uint32_t>(Field::getMaxModulus());
		}

		static integer maxModulus()
		{
			return static_cast<integer>(Field::getMaxModulus());
		}

		static bool goodModulus( const integer& i )
		{
			integer max;
			maxModulus( max );
			if( max == -1 )
				return ( i >= 2 );
			else if( max == 0 )
				return ( i == 0 );
			else
				return ( i >= 2 && i <= max );
		}

		static integer& maxExponent( integer& i )
		{
			return i = 1;
		}

		static bool goodExponent( const integer& i )
		{
			integer max;
			maxExponent( max );
			if( max == -1 )
				return ( i >= 1 );
			else
				return ( i >= 1 && i <= max );
		}
	};


} // Namespace LinBox

namespace LinBox { /*  areFieldEqual  */

	template<class _Field1, class _Field2>
	bool areFieldEqual (const _Field1 &F, const _Field2 &G)
	{
		return false ;
	}

	template<class _Field, class _Category>
	bool areFieldEqualSpecialised(const _Field &F, const _Field &G,
				      const _Category & m)
	{
		return true ;
	}

	template<class _Field>
	bool areFieldEqualSpecialised(const _Field &F, const _Field &G,
				      const RingCategories::ModularTag & m)
	{
		return ( F.characteristic() == G.characteristic() ) ;
	}

	template<class _Field>
	bool areFieldEqual (const _Field &F, const _Field &G)
	{
		return areFieldEqualSpecialised( F,G,typename FieldTraits<_Field>::categoryTag() ) ;
	}

} // Namespace LinBox

#endif // __LINBOX_field_traits_H


