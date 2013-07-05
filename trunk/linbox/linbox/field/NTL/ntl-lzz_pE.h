/* linbox/field/ntl-z_pE.h
 * Copyright (C) 2004  Pascal Giorgi
 * Copyright (C) 2011 LinBox
 *
 * Written by  Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * Modified by W. J. Turner <wjturner@acm.org>
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

/*! @file field/NTL/ntl-lzz_pE.h
 * @ingroup field
 * @ingroup NTL
 * @brief NO DOC
 */

#ifndef __LINBOX_field_ntl_lzz_pe_H
#define __LINBOX_field_ntl_lzz_pe_H

#ifndef __LINBOX_HAVE_NTL
#error "you need NTL here"
#endif

#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pE.h>
#include <time.h>
#include "linbox-config.h"
#include "linbox/util/debug.h"

#include "linbox/field/unparametric.h"
#include "linbox/randiter/unparametric.h"
#include "linbox/field/field-traits.h"



#include "linbox/integer.h"

namespace Givaro
{
	template<>
	NTL::zz_pE& Caster(NTL::zz_pE &x, const Integer &y)
	{
		x=NTL::to_zz_pE(static_cast<long>(y));
		return x;
	}
	template<>
	NTL::zz_pE& Caster(NTL::zz_pE &x, const double &y)
	{
		x=NTL::to_zz_pE(static_cast<long>(y));
		return x;
	}

	template<>
	Integer& Caster (Integer& x, const NTL::zz_pE &y) {
		NTL::zz_pX poly = rep(y);
		Integer base = static_cast<Integer>(NTL::zz_p::modulus());
		long i = deg(poly)+1;
		x = 0;
		for( ; i-- ; ) {
			x *= base;
			x +=  NTL::to_long(rep(coeff(poly, i)));
		}
		return x;
	}
} // namespace Givaro



// Namespace in which all LinBox library code resides
namespace LinBox
{

	//! use ZZ_pEBak mechanism too ?
	class NTL_zz_pE_Initialiser {
	public :
		NTL_zz_pE_Initialiser( const Integer & p, const Integer & k) {
			NTL::zz_p::init( (long) p);
			NTL::zz_pX irredPoly = NTL::BuildIrred_zz_pX ((long) k);
			NTL::zz_pE::init(irredPoly);

		}

		// template <class ElementInt>
		// NTL_zz_pE_Initialiser(const ElementInt& d) {
			// NTL::ZZ_p::init (NTL::to_ZZ(d));
		// }

		// NTL_zz_pE_Initialiser (const NTL::ZZ& d) {
			// NTL::ZZ_p::init(d);
		// }

	};





	/*! @brief zz_pE
	 * Define a parameterized class to easily handle UnparametricField<NTL::zz_pE> field
	 */

	/// \brief for large cardinality, small prime.  \ingroup field
	class NTL_zz_pE : public NTL_zz_pE_Initialiser, public FFPACK::UnparametricOperations<NTL::zz_pE> {
	public:
		typedef NTL::zz_pE Element ;
		typedef FFPACK::UnparametricOperations<Element> Father_t ;
		typedef UnparametricRandIter<Element> RandIter;

		const Element zero,one,mOne ;

		NTL_zz_pE (const integer &p, const integer &k) :
			NTL_zz_pE_Initialiser(p,k),Father_t ()
			,zero( NTL::to_zz_pE(0)),one( NTL::to_zz_pE(1)),mOne(-one)

		{

		}

		Element& random (Element& x) const
		{
			NTL::random(x);
			return x;
		}


		bool isZero (const Element& a) const
		{
			return NTL::IsZero(a);
		}


		bool isOne (const Element& a) const
		{
			return NTL::IsOne(a);
		}

	bool isMOne (const Element& x) const
		{
			Element y ; neg(y,x);
			return isOne(y);
		}


		integer& characteristic (integer &c) const
		{
			return c = static_cast<integer>(NTL::zz_p::modulus());
		}


		integer& cardinality(integer& c) const
		{
			NTL::ZZ card = Element::cardinality();
			long b = NumBytes(card);
			unsigned char* byteArray;
			byteArray = new unsigned char[(size_t)b ];
			BytesFromZZ(byteArray, card, b);

			integer base(256);
			c= integer(0);

			for(long i = b - 1; i >= 0; --i) {
				c *= base;
				c += integer(byteArray[i]);
			}
			delete [] byteArray;

			return c;
		}



		Element& inv(Element& x, const Element& y) const
		{
			x=one/y;
			return x;
		}

		Element& invin(Element& x) const
		{
			x=one/x;
			return x;
		}



		std::istream& read(std::istream& is, Element& x) const
		{
			long tmp = 0;
			is>>tmp;
			x=NTL::to_zz_pE(tmp);
			return is;
		}
	}; // end o class NTL_zz_pE






	template <class Ring>
	struct ClassifyRing;

	template<>
	struct ClassifyRing<UnparametricRandIter<NTL::zz_pE> > {
		typedef RingCategories::ModularTag categoryTag;
	};

	template<>
	class UnparametricRandIter<NTL::zz_pE> {
	public:
		typedef NTL::zz_pE Element;
		UnparametricRandIter<NTL::zz_pE>(const NTL_zz_pE & F ,
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

		UnparametricRandIter<NTL::zz_pE>(const UnparametricRandIter<NTL::zz_pE>& R) :
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

#endif //__LINBOX_ntl_lzz_pe_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

