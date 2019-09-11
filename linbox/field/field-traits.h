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
#include <givaro/givcaster.h>
#include <givaro/givinteger.h>
#include <givaro/givrational.h>
#include <givaro/zring.h>
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include <givaro/gfq.h>

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
		struct GenericTag {
			static std::string name() { return "RingCategories::GenericTag"; }
		};
		//! If it is isomorphic to Z/mZ, for some m or its extensions.
		struct ModularTag : public virtual GenericTag {
			static std::string name() { return "RingCategories::ModularTag"; }
		};
		//! Galois Field  GF(p^e)
		struct GaloisTag : public virtual GenericTag {
			static std::string name() { return "RingCategories::GaloisTag"; }
		};
		//! If it is isomorphic to Z
		struct IntegerTag : public virtual GenericTag {
			static std::string name() { return "RingCategories::IntegerTag"; }
		};
		//! If it is isomorphic to Q
		struct RationalTag : public virtual GenericTag {
			static std::string name() { return "RingCategories::RationalTag"; }
		};
	}

	/*! Default ring category.
	 */

	template <class Field>
	struct ClassifyRing {
		typedef	RingCategories::GenericTag categoryTag;
	};

    template<>
    struct ClassifyRing<Givaro::QField<Givaro::Rational>> {
            typedef RingCategories::RationalTag categoryTag;
    };
    template<>
    struct ClassifyRing<Givaro::ZRing<Givaro::Integer>> {
            typedef RingCategories::IntegerTag categoryTag;
    };
//     template<>
//     struct ClassifyRing<Givaro::IntegerDom> {
//             typedef RingCategories::IntegerTag categoryTag;
//     };


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

	/*! FieldTrait.
	 */
	template <class Field>
	struct FieldTraits {

		typedef typename ClassifyRing<Field>::categoryTag categoryTag;

        typedef std::is_same<categoryTag,RingCategories::ModularTag> is_modular;

		static integer& maxModulus( integer& i )
		{
            		return Caster(i, Field::maxCardinality());
		}

		static uint64_t & maxModulus( uint64_t& i )
		{
            		return Caster(i, Field::maxCardinality());
// 			return i = static_cast<uint64_t>(Field::maxCardinality());
		}

		static uint32_t & maxModulus( uint32_t& i )
		{
            		return Caster(i, Field::maxCardinality());
// 			return i = static_cast<uint32_t>(Field::maxCardinality());
		}

		template<class T>
		static T & maxModulus( T & i )
		{
			return i = static_cast<T>(Field::maxCardinality());
		}

		static integer maxModulus()
		{
            		integer maxm;
            		return Caster(maxm, Field::maxCardinality());
// 			return static_cast<integer>(Field::maxCardinality());
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

                    /* \brief returns the best modulus bitsize to use for e.g. ChineseRemaindering
                     */
                static inline uint64_t bestBitSize(){return maxModulus().bitsize()-1;}

                    /* \brief returns the best modulus bitsize to use for e.g. ChineseRemaindering,
                     * given the dimension of linear algebra operations to be performed.
                     * Will be specialized for fields with delayed modular reductions.
                     */
                static inline uint64_t bestBitSize(size_t n){return bestBitSize();}
	};

        template<>
        inline uint64_t FieldTraits<Givaro::Modular<double> >::bestBitSize(size_t n){return std::max (UINT64_C(22), uint64_t(52-(n?log2(n):0))>>1); }
        template<>
        inline uint64_t FieldTraits<Givaro::ModularBalanced<double> >::bestBitSize(size_t n){return std::max (UINT64_C(23), uint64_t(54-(n?log2(n):0))>>1);}
        template<>
        inline uint64_t FieldTraits<Givaro::Modular<float> >::bestBitSize(size_t n){return std::max (UINT64_C(8), uint64_t(24-(n?log2(n):0))>>1);}
        template<>
        inline uint64_t FieldTraits<Givaro::ModularBalanced<float> >::bestBitSize(size_t n){return std::max ( UINT64_C(8), uint64_t(26-(n?log2(n):0))>>1);}
        template<>
        inline uint64_t FieldTraits<Givaro::Modular<int64_t> >::bestBitSize(size_t n){return std::max ( UINT64_C(26), uint64_t(63-(n?log2(n):0))>>1);}
        template<>
        inline uint64_t FieldTraits<Givaro::ModularBalanced<int64_t> >::bestBitSize(size_t n){return std::max ( UINT64_C(27), uint64_t(64-(n?log2(n):0))>>1);}

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

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
