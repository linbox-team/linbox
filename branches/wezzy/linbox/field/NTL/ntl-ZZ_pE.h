/* Copyright (C) LinBox
 *
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

/*! @file field/NTL/ntl-ZZ_pE.h
 * @ingroup field
 * @ingroup NTL
 * @brief NO DOC
 */

#ifndef __LINBOX_field_ntl_zz_pe_H
#define __LINBOX_field_ntl_zz_pe_H

#ifndef __LINBOX_HAVE_NTL
#error "you need NTL here"
#endif

#include <time.h>
#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ.h>

#include "linbox/field/unparametric.h"
#include "linbox/randiter/unparametric.h"
#include "linbox/field/field-traits.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#endif //__LINBOX_XMLENABLED



#include <givaro/givcaster.h>
namespace Givaro
{

	template<>
	NTL::ZZ_pE& Caster(NTL::ZZ_pE &x, const Integer &y)
	{
		x=NTL::to_ZZ_pE(static_cast<long>(y));
		return x;
	}
	template<>
	NTL::ZZ_pE& Caster(NTL::ZZ_pE &x, const double &y)
	{
		x=NTL::to_ZZ_pE(static_cast<long>(y));
		return x;
	}

	// Rich Seagraves, 7-15-03
	// On the orders of Dr Saunders, I'm re-writing init & convert so that
	// they convert a ZZpE into a padic number, ie a0 + a1x + a2x^2 +... ->
	// a0 + a1*p + a2*p^2 + ...
	//

	Integer& Caster(Integer& c, const NTL::ZZ_pE& e)
	{
		NTL::ZZ_pX poly = rep(e);
		//Integer base = _p;
		Integer base = static_cast<Integer>(to_long(NTL::ZZ_p::modulus()));
		long i;

		c = 0;
		for(i = deg(poly); i >= 0; --i) {
			c *= base;
			c +=  NTL::to_long(rep(coeff(poly, i)));
		}

		return c;
	}
} // namespace Givaro

// Namespace in which all LinBox library code resides
namespace LinBox
{

	class NTL_ZZ_pE_Initialiser {
	public :
		NTL_ZZ_pE_Initialiser( const Integer & p, size_t k ) {
			// linbox_check(e == 1);
			// if ( q > 0 )
			// NTL::ZZ_p::init(NTL::to_ZZ((std::string(q)).data())); // it's an error if q not prime, e not 1
			NTL::ZZ_p::init(NTL::to_ZZ(std::string(p).data()));
			NTL::ZZ_pX irredPoly = NTL::BuildIrred_ZZ_pX ((long) k);
			NTL::ZZ_pE::init(irredPoly);

		}

		// template <class ElementInt>
		// NTL_ZZ_pE_Initialiser(const ElementInt& d) {
			// NTL::ZZ_p::init (NTL::to_ZZ(d));
		// }

		// NTL_ZZ_pE_Initialiser (const NTL::ZZ& d) {
			// NTL::ZZ_p::init(d);
		// }

	};



	/**
	 * @brief Wrapper of ZZ_pE from NTL
	 * Define a parameterized class to handle easily UnparametricField<NTL::ZZ_pE> field
	 */
	class NTL_ZZ_pE : public NTL_ZZ_pE_Initialiser, public UnparametricOperations<NTL::ZZ_pE> {
	public:
		typedef NTL::ZZ_pE Element ;
		typedef FFPACK::UnparametricOperations<Element> Father_t ;

		const Element zero,one,mOne ;


		NTL_ZZ_pE (const integer &p, const integer &k) :
			NTL_ZZ_pE_Initialiser(p,k),Father_t ()
			,zero( NTL::to_ZZ_pE(0)),one( NTL::to_ZZ_pE(1)),mOne(-one)


		{

				}


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
		return c=static_cast<integer>(to_long(NTL::ZZ_p::modulus()));
		//NTL::ZZ_p::modulus();
	}


	integer& cardinality(integer& c) const
	{
		c=static_cast<integer>(to_long(NTL::ZZ_p::modulus()));
		c=pow(c,Element::degree());
		return c;
	}

	Element& inv(Element& x, const Element& y) const
	{
		x=NTL::to_ZZ_pE(1)/y;
		return x;
	}

	Element& invin(Element& x) const
	{
		x=NTL::to_ZZ_pE(1)/x;
		return x;
	}


	std::istream& read(std::istream& is, Element& x) const
	{
		long tmp;
		is>>tmp;
		x=NTL::to_ZZ_pE(tmp);
		return is;
	}



#ifdef __LINBOX_XMLENABLED


	bool toTag(LinBox::Writer &W) const
	{
		std::string s;
		NTL::ZZ_pX poly = Element::modulus();
		long i;

		W.setTagName("field");
		W.setAttribute("implDetail", "ntl-ZZpE");
		W.setAttribute("cardinality", LinBox::Writer::numToString(s, _card));

		W.addTagChild();
		W.setTagName("finite");

		W.addTagChild();
		W.setTagName("characteristic");
		W.addNum(_p);
		W.upToParent();

		W.addTagChild();
		W.setTagName("extension");
		W.addNum(deg(poly) + 1);
		W.upToParent();

		W.addTagChild();
		W.setTagName("polynomial");

		for(i = 0; i <= deg(poly); ++i) {
			W.addNum(coeff(poly, i));
		}
		W.upToParent();
		W.upToParent();
		W.upToParent();

		return true;
	}


	std::ostream &write(std::ostream &os) const
	{
		LinBox::Writer W;
		if( toTag(W) )
			W.write(os);

		return os;
	}


	// Elemnt Reading & writing functions
	// BIG NOTE:  It was decided that for extension fields, the elements
	// would be represented using a single number that has the following
	// property:  for an element e in ZZp[x], with e = a0 + a1x + a2x^2 + ...,
	// represent e as "<cn>n</cn>" where n = a0 + a1 * p + a2 * p^2 + ...
	//


	bool toTag(LinBox::Writer &W, const Element &e) const
	{
		NTL::ZZ_pX poly = rep(e);
		NTL::ZZ accum, base = NTL::ZZ_p::modulus();
		long i;
		std::string s;

		accum = 0;
		for(i = deg(poly); i >= 0; --i) {
			accum *= base;
			accum += rep(coeff(poly, i));
		}


		W.setTagName("cn");
		W.addDataChild(LinBox::Writer::numToString(s, accum));

		return true;
	}


	std::ostream &write(std::ostream &os, const Element &e) const
	{

		LinBox::Writer W;
		if( toTag(W, e))
			W.write(os);

		return os;
	}




	bool fromTag(LinBox::Reader &R, Element &e) const
	{
		NTL::ZZ total, base = NTL::ZZ_p::modulus(), rem;
		std::stringstream ss;

		if(!R.expectTagName("cn") || !R.expectChildTextNum(total))
			return false;

		ss << "[";
		while(total > 0) {
			rem = total % base;
			total /= base;
			ss << rem;
			if(total > 0) ss << " ";
		}

		ss << "]";

		ss >> e; // use the extraction stream operator

		return true;
	}


	std::istream &read(std::istream &is, Element &e) const
	{
		LinBox::Reader R(is);
		if( !fromTag(R, e)) {
			is.setstate(std::istream::failbit);
			if(!R.initalized()) {
				is.setstate(std::istream::badbit);
			}
		}

		return is;
	}


#endif


	}; // end of class NTL_ZZ_pE






	template <class Ring>
	struct ClassifyRing;

	template<>
	struct ClassifyRing<UnparametricRandIter<NTL::ZZ_pE> > {
		typedef RingCategories::ModularTag categoryTag;
	};



}

namespace LinBox
{


	template<>
	class UnparametricRandIter<NTL::ZZ_pE> {
	public:
		typedef NTL::ZZ_pE Element;
		UnparametricRandIter<NTL::ZZ_pE>(const NTL_ZZ_pE & F ,
						 const integer& size =0,
						 const integer& seed =0
						) :
			_size(size), _seed(seed)
		{
			if(_seed == 0)
				NTL::SetSeed(NTL::to_ZZ(time(0)));
			else
				NTL::SetSeed( NTL::to_ZZ(static_cast<long>(seed)) );
		}

#ifdef __LINBOX_XMLENABLED
		UnparametricRandIter<NTL::ZZ_pE>(LinBox::Reader &R)
		{
			if(!R.expectTagName("randiter")) return;
			if(!R.expectAttributeNum("seed", _seed) || !R.expectAttributeNum("size", _size)) return;

			if(_seed == 0) _seed = time(NULL);

			NTL::SetSeed(NTL::to_ZZ(_seed));
		}
#endif


		UnparametricRandIter<NTL::ZZ_pE>(const UnparametricRandIter<NTL::ZZ_pE>& R) :
			_size(R._size), _seed(R._seed)

		{
			if(_seed == 0)
				NTL::SetSeed(NTL::to_ZZ(time(0)));
			else
				NTL::SetSeed(NTL::to_ZZ( static_cast<long>(_seed)) );
		}

		NTL::ZZ_pE& random (NTL::ZZ_pE& x) const
		{
			NTL::random(x);
			return x;
		}

	protected:
		size_t _size;
		size_t _seed;
	};
}
#endif //__LINBOX_field_ntl_zz_pe_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

