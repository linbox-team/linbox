/* Copyright (C) 2010 LinBox
 * written by bds, wan
 *
 *
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
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef __LINBOX_local2_32_H
#define __LINBOX_local2_32_H

#include "linbox/field/unparametric.h"
#include "linbox/util/debug.h"
#include "linbox/linbox-config.h"
#include <linbox/field/field-traits.h>
#include <linbox/integer.h>
#include <linbox/field/field-traits.h>

namespace LinBox
{
  
        template<typename Ring>
	struct ClassifyRing;

	class Local2_32;

	template<>
	struct ClassifyRing<Local2_32> {
		typedef RingCategories::ModularTag categoryTag;
	};

	/** \brief Fast arithmetic mod 2^32, including gcd.
	 *
	 * Extend UnparametricField<uint32> which is a representation 
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

	struct Local2_32: public UnparametricField<uint32>
	{
	public:

		typedef UnparametricField<uint32>::Element Element;
		typedef enum {_min=0,_max=32} Exponent; // enum?
		//Exponent& init(Exponent& a) { return a = 32; }
			
		Local2_32 (int p=2, int exp=32) :UnparametricField<uint32>(p,exp) {
			if(p != 2) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be 2");
			if(exp != 32) throw PreconditionFailed(__FUNCTION__,__LINE__,"exponent must be 32");
		}

		/*
		static inline Element& gcd(Element& c, Element& a, const Element& b)
		{   c = a | b; Exponent k = 0; 
		    while (! (c & 1)) {c >>= 1; ++k;}
		    //gcdin (k, b);
		    cout << "gcd called" << endl;
		    return c = 1 << k;
		}

		*/
		// assume k is an exponent of 2.
		static inline Exponent& gcdin(Exponent& k, const Element& b)
		{   /*
		      Element c = b >> k;
		      c <<= k;
		      Element d = b;
		      std::cout << "c, b" << c << " " << b <<  "\n";
		      if (c != b)  for(k = 0; ! (d & 1); ++k) d >>= 1; 
		      std::cout << "g, b =" << (int)k << " " << b << "\n";
		    */
			Element d = b;
			int i;
			for ( i = 0; (i < k) && (!(d&1)); ++ i) d >>= 1;
			return k = Exponent(i);
		}

		static inline bool isUnit(const Exponent& a)
		{   return a == 0;   }

		static inline bool isZero(const Element& a)
		{   return a == 0;   }

		static inline bool isZero(const Exponent& a)
		{   return a >= 32;   }

		// not used ...
		static inline bool isUnit(const Element& a)
		{   return a & 1;   }

		//Element& div(Element& c, const Element& a, const Element& b) const
		//{   return c = NTL::rep(a)/NTL::GCD(NTL::rep(a),NTL::rep(b));   }
		//
		
		static inline Element& mulin(Element& a, const Exponent& k) 
		{  
			if (k >= 32) return a = 0;
			else return a <<= k;  
		}

		static inline Element& mulin(Element& a, const Element& b)  {
			return a *= b;
		}

		static inline Element& axpyin(Element& r, const Element& x, const Element& y) {
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
		static inline Element& divin(Element& a, const Exponent& k)
		{   return a >>= k;   }

		static inline Element& inv(Element& a, const Element& b) {

			if (!isUnit(b))
				throw PreconditionFailed(__FUNCTION__,__LINE__,"inv: not a unit");
			else {

				Element g, s, t;
				
				xgcd(g, s, t, b, -b);
				
				return a = s - t;
			}
		}

                static inline integer getMaxModulus()
                        { return integer( "4294967296" ); } // 2^32

	protected:
		
		static  Element& xgcd(Element& d, Element& s, Element& t, const Element& a, const Element& b) 
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
			//std::cout << "XGCD is called: d, s, t, a, b, sa + tb: " << d << ' '
			//	<< s << ' ' << t << ' ' << a << ' ' << b << ' ' << s * a + t * b << '\n';
			return d;

			/*

			//Element  u, v, u0, v0, u1, v1, u2, v2, q, r;
			
			Element u, v, q, r;
			
			int64  u0, u1, u2;
			
			u1 = 1; //v1 = 0;
			u2 = 0; //v2 = 1;
			u = a; v = b;
		    
			if ( b == 0) {
				s = 1;
				t = 0;
				return d = a ;
			}
			
			if (v != 0) {
				q = u / v;
				//r = u % v;
				r = u - q*v;
				u = v;
				v = r;
				u0 = u2;
				//v0 = v2;
				u2 =  u1 - q * u2;
				//v2 = v1- q * v2;
				u1 = u0;
				//v1 = v0;
				
			}

			while (v != 0) {
				r = u;
				while ( r >= v) {
					r = u - v;
					
					u2 = u1 - u2;
				}
				u0 = u2;
				u1 = u0;
				u = v;
				v = r;
				
				
			}	
				
			
			while (v != 0) {
				q = u / v;
				//r = u % v;
				r = u - q*v;
				u = v;
				v = r;
				
				u0 = u2;
				//v0 = v2;
				u2 =  u1 - q * u2;
				//v2 = v1- q * v2;
				u1 = u0;
				//v1 = v0;
			}
			
		    
			d = u;
			s = u1;
			
			t = ((int64) d - u1 * (int64) a) / (int64)b;

			//std::cout << "XGCD is called: d, s, t, a, b, sa + tb: " << d << ' '
			//<< s << ' ' << t << ' ' << a << ' ' << b << ' ' << s * a + t * b << '\n';
			return d;
			
			
			*/
		}
		

		/** @brief 
		 * Half GCD
		 * g = gcd (a, b).
		 * exists t, such that: s * a + t * b = g.
		 * return g.
		 */
		static Element& HGCD (Element& g, Element& s, const Element& a, const  Element& b) {
			
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
		return i == Local2_32::getMaxModulus();
	}
		
} // namespace LinBox

#endif // __LINBOX_local2_32_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
