/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/field/ntl-z_pE.h
 * Copyright (C) 2004  Pascal Giorgi
 * Copyright (C) 2011 LinBox
 *
 * Written by  Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * Modified by W. J. Turner <wjturner@acm.org>
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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/*! @file field/ntl-GF2E.h
 * @ingroup field
 * @ingroup NTL
 * @brief NO DOC
 */

#ifndef __LINBOX_field_ntl_gf2e_H
#define __LINBOX_field_ntl_gf2e_H

#ifndef __LINBOX_HAVE_NTL
#error "you need NTL here"
#endif

#include <NTL/GF2XFactoring.h>
#include <NTL/GF2E.h>
#include <time.h>

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#include "linbox/field/unparametric.h"
#include "linbox/randiter/unparametric.h"
#include "linbox/field/field-traits.h"


// Namespace in which all LinBox library code resides
namespace LinBox
{
	template<>
	//	NTL::GF2E& UnparametricField<NTL::GF2E>::init (NTL::GF2E &x, const integer &y) const
	NTL::GF2E& Caster(NTL::GF2E &x, const integer &y)
	{
		x=NTL::to_GF2E(static_cast<long>(y));
		return x;
	}
	template<>
	//	NTL::GF2E& UnparametricField<NTL::GF2E>::init (NTL::GF2E &x, const double &y) const
	NTL::GF2E& Caster(NTL::GF2E &x, const double &y)
	{
		x=NTL::to_GF2E(static_cast<long>(y));
		return x;
	}


	template<>
	//	integer& UnparametricField<NTL::GF2E>::convert (integer& x, const NTL::GF2E &y) const	{
	integer& Caster(integer& x, const NTL::GF2E &y)
	{
		NTL::GF2X poly = rep(y);

		long i;
		x = 0;
		for(i = deg(poly); i >= 0; --i) {
			x <<= 1;
			x += rep(coeff(poly, i));
		}
		return x;
	}



	/// \ingroup field

	class NTL_GF2E_Initialiser {
	public :
		NTL_GF2E_Initialiser( const Integer & p, size_t k = 1) {
			linbox_check(p == 2);
			// if(p != 2) throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus must be 2");
			NTL::GF2X irredPoly = NTL::BuildSparseIrred_GF2X((long) k);
			NTL::GF2E::init(irredPoly);

		}

		// template <class ElementInt>
		// NTL_GF2E_Initialiser(const ElementInt& d) {
			// NTL::ZZ_p::init (NTL::to_ZZ(d));
		// }

		// NTL_GF2E_Initialiser (const NTL::ZZ& d) {
			// NTL::ZZ_p::init(d);
		// }

	};


	/*
	 * Define a parameterized class to easily handle UnparametricField<NTL::GF2E> field
	 */

	class NTL_GF2E :  public NTL_GF2E_Initialiser, public FFPACK::UnparametricOperations<NTL::GF2E> {
	public:
		typedef NTL::GF2E Element ;
		typedef FFPACK::UnparametricOperations<Element> Father_t ;
		typedef UnparametricRandIter<Element> RandIter;

		const Element zero,one,mone ;

		NTL_GF2E (const integer &p, const integer &k) :
			NTL_GF2E_Initialiser(p,k),Father_t ()
			,zero( NTL::to_GF2E(0)),one( NTL::to_GF2E(1)),mone(-one)
		{ }

		bool isZero (const Element& a) const
		{
			return NTL::IsZero(a);
		}


		bool isOne (const Element& a) const
		{
			return NTL::IsOne(a);
		}

		integer& characteristic (integer &c) const
		{
			return c = 2;
		}

		integer& cardinality(integer& c) const
		{
			c=1;
			c<<= Element::degree();
			return c;
		}

		Element& inv(Element& x, const Element& y) const
		{
			x=NTL::to_GF2E(1)/y;
			return x;
		}

		Element& invin(Element& x) const
		{
			x=NTL::to_GF2E(1)/x;
			return x;
		}

		std::istream& read(std::istream& is, Element& x) const
		{
			long tmp;
			is>>tmp;
			x=NTL::to_GF2E(tmp);
			return is;
		}
	}; // end o class NTL_GF2E

	template <class Ring>
	struct ClassifyRing;

	template <class Element>
	struct ClassifyRing<UnparametricRandIter<Element> >;

	template <>
	struct ClassifyRing<UnparametricRandIter<NTL::GF2E> > {
		typedef RingCategories::ModularTag categoryTag;
	};

	template<>
	class UnparametricRandIter<NTL::GF2E> {
	public:
		typedef NTL::GF2E Element;
		UnparametricRandIter<NTL::GF2E>(const NTL_GF2E & F,
						const size_t& size = 0,
						const size_t& seed = 0
					       ) :
			_size(size), _seed(seed)
		{
			if(_seed == 0)
				NTL::SetSeed(NTL::to_ZZ(time(0)));
			else
				NTL::SetSeed(NTL::to_ZZ(_seed));
		}

		UnparametricRandIter<NTL::GF2E>(const UnparametricRandIter<NTL::GF2E>& R) :
			_size(R._size), _seed(R._seed)

		{
			if(_seed == 0)
				NTL::SetSeed(NTL::to_ZZ(time(0)));
			else
				NTL::SetSeed(NTL::to_ZZ(_seed));
		}

		Element& random (Element& x) const
		{
			NTL::random(x);
			return x;
		}

	protected:
		size_t _size;
		size_t _seed;
	};


}

#endif //__LINBOX_ntl_gf2e_H

