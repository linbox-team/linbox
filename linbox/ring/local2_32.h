/* Copyright (C) 2010 LinBox
 * written by bds, wan
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


#ifndef __LINBOX_local2_32_H
#define __LINBOX_local2_32_H

#include <givaro/zring.h>
#include "linbox/util/debug.h"
#include "linbox/linbox-config.h"
#include "linbox/field/field-traits.h"
#include "linbox/integer.h"
#include "linbox/field/hom.h"

namespace LinBox
{

	template<typename Ring>
	struct ClassifyRing;

	struct Local2_32;

	template<>
	struct ClassifyRing<Local2_32> {
		typedef RingCategories::ModularTag categoryTag;
	};

	/** \brief Fast arithmetic mod 2^32, including gcd.
	 *
	 * Extend Givaro::ZRing<uint32_t> which is a representation
	 * of Z_2^32. It is especially fast because it uses hardware arithmetic
	 * directly.  This ring is a Local Principal Ideal Ring.
	 *
	 * These needed PIR functions are added:
	 * gcdin(), isUnit(), also inv() is modified to work correctly.
	 * The type Exponent is added: more effective rep of the powers of 2,
	 * which are important because gcds are powers of 2).
	 * This entails some new versions of divin(), mulin(), isUnit().
	 *
	 * Those are the function needed for the LocalSmith algorithm.
	 * Further appropriate PIR functions may be added later.
	 * \ingroup field
	 */

	struct Local2_32: public Givaro::ZRing<uint32_t>
	{
	public:

		typedef Givaro::ZRing<uint32_t>::Element Element;
		typedef enum {_min=0,_max=32} Exponent; // enum?
		//Exponent& init(Exponent& a) { return a = 32; }

		Local2_32 (int p=2, int exp=32) :
			Givaro::ZRing<uint32_t>()
		{
			if(p != 2) throw PreconditionFailed(LB_FILE_LOC,"modulus must be 2");
			if(exp != 32) throw PreconditionFailed(LB_FILE_LOC,"exponent must be 32");
		}

		inline Element& gcd(Element& c, const Element& a, const Element& b) const
		{   c = a | b; 
			if (c == 0) return c;
			uint32_t i = 0;
		    while (! (c & 1)) {c >>= 1; ++i;}
			return c = 1 << i;
		}

		inline Element& gcdin(Element& b, const Element& a) const
		{	
			Element c = b; return gcd(b, c, a); }
			/*
			if (isZero(b)) return b = a;
			Element d = b;
			Exponent k;
			int32_t i;
			for ( i = 0; (i < k) && (!(d&1)); ++ i) d >>= 1;
			k = Exponent(i);
			gcdin(k, a);
			return b = 1 << k;
		}
		*/
		// assume k is an exponent of 2.
		inline Exponent& gcdin(Exponent& k, const Element& b) const
		{   
			Element d = b;
			int32_t i;
			for ( i = 0; (i < k) && (!(d&1)); ++ i) d >>= 1;
			return k = Exponent(i);
		}

		inline bool isUnit(const Element& a) const
		{   return a & 1;   }

		inline bool isUnit(const Exponent& a) const
		{   return a == 0;   }

		inline bool isZero(const Element& a) const
		{   return a == 0;   }

		inline bool isZero(const Exponent& a) const
		{   return a >= 32;   }

		//Element& div(Element& c, const Element& a, const Element& b) const
		//{   return c = NTL::rep(a)/NTL::GCD(NTL::rep(a),NTL::rep(b));   }
		//

		inline Element& mulin(Element& a, const Exponent& k) const
		{
			if (k >= 32) return a = 0;
			else return a <<= k;
		}

		inline Element& mulin(Element& a, const Element& b)  const {
			return a *= b;
		}

		inline Element& axpyin(Element& r, const Element& x, const Element& y)  const{
			return r += x * y;
		}

		/*
		   static inline bool isDivisor(Element a, Element b)
		   {   while (! (a ^ 1))
		   {   if (b ^ 1) return false;
		   a = a >> 1; b = b >> 1;
		   }
		   return true;
		   }
		   */

		// assume k is an exponent of 2 and the power of 2 exactly divides a
		inline Element& divin(Element& a, const Exponent& k) const
		{   return a >>= k;   }

		inline Element& inv(Element& a, const Element& b) const {

			if (!isUnit(b))
				throw PreconditionFailed(LB_FILE_LOC,"inv: not a unit");
			else {

				Element g, s, t;

				xgcd(g, s, t, b, -b);

				return a = s - t;
			}
		}

		static inline integer maxCardinality()
		{ return integer( "4294967296" ); } // 2^32

	protected:

		Element& xgcd(Element& d, Element& s, Element& t, const Element& a, const Element& b) const
		{

			Element  u, v, u0, v0, u1, v1, u2, v2, q, r;

			u1 = 1; v1 = 0;
			u2 = 0; v2 = 1;
			u = a; v = b;

			while (v != 0) {
				q = u / v;
				//r = u % v;
				r = u - q*v;
				u = v;
				v = r;
				u0 = u2;
				v0 = v2;
				u2 =  u1 - q*u2;
				v2 = v1- q*v2;
				u1 = u0;
				v1 = v0;
			}


			d = u;
			s = u1;
			t = v1;
			return d;

		}


		/** @brief
		 * Half GCD
		 * g = gcd (a, b).
		 * exists t, such that: s * a + t * b = g.
		 * return g.
		 */
		Element& HGCD (Element& g, Element& s, const Element& a, const  Element& b) const {

			Element  u, v, u0, u1, u2, q, r;

			u1 = 1;
			u2 = 0;
			u = a; v = b;

			while (v != 0) {
				q = u / v;
				//r = u % v;
				r = u - q*v;
				u = v;
				v = r;

				u0 = u2;

				u2 =  u1 - q*u2;

				u1 = u0;

			}


			g = u;
			s = u1;

			return g;

		}

	};

	template<>
	bool FieldTraits< Local2_32 >::goodModulus( const integer& i ) {
		return i == Local2_32::maxCardinality();
	}

	template<>
	class Hom<Givaro::ZRing<Integer>, Local2_32> {
	public:
		typedef Givaro::ZRing<Integer> Source;
		typedef Local2_32 Target;
		typedef typename Source::Element SrcElt;
		typedef typename Target::Element Elt;

		Hom(const Source& S, const Target& T) :
			_source (S), _target (T)
		{}
		inline Elt& image(Elt& t, const SrcElt& s)
		{
			SrcElt x = s % Local2_32::maxCardinality();
			if (x < 0) x += Local2_32::maxCardinality();
			_target. init (t, x);
			return t;
		}
		inline SrcElt& preimage(SrcElt& s, const Elt& t)
		{
			_target. convert (s, t);
			return s;
		}
		const Source& source() { return _source;}
		const Target& target() { return _target;}

	protected:
		const Source& _source;
		const Target& _target;
	}; // end Hom

} // namespace LinBox

#endif // __LINBOX_local2_32_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
