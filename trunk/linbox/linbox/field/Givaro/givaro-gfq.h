/* linbox/field/givaro-gfq.h
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 * JGD 12.06.2002 : -- I don't see the need of *(new in convert
 * JGD 19.09.2003 : added isZero
 * WJT 24.06.2005 : Removed using declarations
 *
 * ------------------------------------
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file field/Givaro/givaro-gfq.h
 *
 * @warning This wrapper works only with an improved version of Givaro ;
 * This version of givaro won't be available for public yet ;
 * But it is available on my web page ;
 * You can send me a mail to get it or for others details.
 */

#ifndef __LINBOX_field_givaro_gfq_H
#define __LINBOX_field_givaro_gfq_H


#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/integer.h"
#include "linbox/field/field-traits.h"
#include "linbox/field/field-interface.h"


//------------------------------------
// Files of Givaro library


#include <givaro/givtablelimits.h>
#include <givaro/givgfq.h>
#include <givaro/giv_randiter.h>
#include <givaro/givpoly1factor.h>
//------------------------------------

// #include "linbox/vector/blas-vector.h"

// Namespace in which all LinBox code resides
namespace LinBox
{
	template<class A, class B>
		class BlasVector ;

	template <class Ring>
	struct ClassifyRing;

	struct PID_integer;

	class GivaroGfq;

	template<>
	struct ClassifyRing<GivaroGfq> {
		typedef RingCategories::ModularTag categoryTag;
	};


	template<>
	inline integer& FieldTraits<GivaroGfq>::maxModulus( integer& i )
	{
		return i = integer( 32749 );  // prevprime( 2^15 )
	}

	template<>
	inline bool FieldTraits<GivaroGfq>::goodModulus( const integer& i )
	{
		integer max;
		if( i < 2 || i > FieldTraits<GivaroGfq>::maxModulus(max) )
			return false;
		return probab_prime( i, 10 );
	}

	template<>
	inline integer& FieldTraits<GivaroGfq>::maxExponent( integer& i )
	{
		return i = _GIVARO_FF_MAXEXPONENT_;  // Cardinality must be <= 2^21
	}


	/** Wrapper of Givaro's GFqDom<int32_t>  class.
	  \ingroup field

	 *  This class allows to construct only extension fields with a prime characteristic.
	 */
	class GivaroGfq : public Givaro::GFqDom<int32_t>, public FieldInterface {

	public:

		typedef Givaro::GFqDom<int32_t> Father_t ;

		using Father_t::one ;
		using Father_t::zero;
		using Father_t::mOne;
		/** Element type.
		 *  This type is inherited from the Givaro class GFqDom<int32_t>
		 */
		typedef  Givaro::GFqDom<int32_t>::Rep Element;

		/** RandIter type
		 *  This type is inherited from the Givaro class GFqDom<TAG>
		 */
		typedef Givaro::GIV_randIter< Givaro::GFqDom<int32_t>, LinBox::integer >  RandIter;

		/** Empty Constructor
		*/
		GivaroGfq() :
		 Givaro::GFqDom<int32_t>()
		{
		}

		/** Constructor from an integer
		 *  this constructor use the ZpzDom<TAG> constructor
		 */
		GivaroGfq(const integer& p, const integer& k=1) :
		 Givaro::GFqDom<int32_t>(static_cast<UTT>(int32_t(p)), static_cast<UTT>(int32_t(k)))
		{
			//enforce that the cardinality must be <2^16, for givaro-gfq
			int32_t pl=p;
			for(int32_t i=1;i<k;++i) pl*=(int32_t)p;
			if(!FieldTraits<GivaroGfq>::goodModulus(p)) {
				throw PreconditionFailed(LB_FILE_LOC,"modulus be between 2 and 2^15 and prime");
			}
			else if(pl>(1<<_GIVARO_FF_MAXEXPONENT_))  {
                            std::stringstream MSGFFMAXE;
                            MSGFFMAXE << "cardinality must be < 2^" << _GIVARO_FF_MAXEXPONENT_;
				throw PreconditionFailed(LB_FILE_LOC,MSGFFMAXE.str().c_str());
			}

		}

		// This constructor takes a vector of ints that represent the polynomial
		// to use (for modular arithmetic on the extension field).
		GivaroGfq(const integer& p, const integer& k, const std::vector<integer>& modPoly) :
		 Givaro::GFqDom<int32_t>(static_cast<UTT>(int32_t(p)), static_cast<UTT>(int32_t(k)), modPoly)
		{
		}

		GivaroGfq(const integer& p, const integer& k, const BlasVector<PID_integer,std::vector<integer> >& modPoly) :
			Givaro::GFqDom<int32_t>(static_cast<UTT>(int32_t(p)), static_cast<UTT>(int32_t(k)), modPoly)
		{
		}

		/** Characteristic.
		 * Return integer representing characteristic of the domain.
		 * Returns a positive integer to all domains with finite characteristic,
		 * and returns 0 to signify a domain of infinite characteristic.
		 * @return integer representing characteristic of the domain.
		 */
		integer& characteristic(integer& c) const
		{
			return c=integer(static_cast<int32_t>( Givaro::GFqDom<int32_t>::characteristic()));
		}

		int32_t characteristic() const
		{
			return static_cast<int32_t>( Givaro::GFqDom<int32_t>::characteristic());
		}

#if (GIVARO_VERSION<30403)
		unsigned long characteristic(unsigned long & c) const
		{
			return c = static_cast<int32_t>( Givaro::GFqDom<int32_t>::characteristic());
		}
#else
		unsigned long characteristic(unsigned long & c) const
		{
			return Givaro::GFqDom<int32_t>::characteristic(c);
		}
#endif

		/** Cardinality.
		 * Return integer representing cardinality of the domain.
		 * Returns a non-negative integer for all domains with finite
		 * cardinality, and returns -1 to signify a domain of infinite
		 * cardinality.
		 * @return integer representing cardinality of the domain
		 */
		integer& cardinality(integer& c) const
		{
			return c=integer(static_cast<int32_t>( Givaro::GFqDom<int32_t>::size()));
		}


		integer cardinality() const
		{
			return integer(static_cast<int32_t>( Givaro::GFqDom<int32_t>::cardinality()));
		}


		/** Initialization of field base Element from an integer.
		 * Behaves like C++ allocator construct.
		 * This function assumes the output field base Element x has already been
		 * constructed, but that it is not already initialized.
		 * We assume that the type of Element is short int.
		 * this methos is just a simple cast.
		 * @return reference to field base Element.
		 * @param x field base Element to contain output (reference returned).
		 * @param y integer.
		 */
		Element& init(Element& x , const integer& y = 0) const
		{
			return Givaro::GFqDom<int32_t>::init( x, int32_t(y % (integer) _q));
		}

		// TO BE OPTIMIZED
		Element& init(Element& x , const float y) const
		{
			return Givaro::GFqDom<int32_t>::init( x, (double)y);
		}

		template<class YYY>
		Element& init(Element& x , const YYY& y) const
		{
			return Givaro::GFqDom<int32_t>::init( x, y);
		}

		/** Conversion of field base Element to an integer.
		 * This function assumes the output field base Element x has already been
		 * constructed, but that it is not already initialized.
		 * @return reference to an integer.
		 * @param x integer to contain output (reference returned).
		 * @param y constant field base Element.
		 */
		integer& convert(integer& x, const Element& y) const
		{
			int32_t tmp;
			return x = integer( Givaro::GFqDom<int32_t>::convert(tmp,y));
		}
		// TO BE OPTIMIZED
		float& convert(float& x, const Element& y) const
		{
			double tmp;
		 Givaro::GFqDom<int32_t>::convert( tmp, y);
			return x = (float)tmp;
		}

		template<class XXX>
		XXX& convert(XXX& x, const Element& y) const
		{
			return Givaro::GFqDom<int32_t>::convert( x, y);
		}

		//bool isZero(const Element& x) const { return Givaro::GFqDom<int32_t>::isZero(x); }


	}; // class GivaroGfq





} // namespace LinBox

#endif // __LINBOX_field_givaro_gfq_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

