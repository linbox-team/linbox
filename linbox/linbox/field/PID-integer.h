/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/field/PID-integer.h
 * Copyright (C) 2004 Pascal Giorgi 
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
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




#ifndef __PID_INTEGER_H
#define __PID_INTEGER_H

#include <linbox/integer.h>
#include <linbox/field/unparametric.h>
#include <linbox/field/field-traits.h>



namespace LinBox {

	template <class Ring>
	class ClassifyRing;

	
	class PID_integer : public UnparametricField<integer> 
	{

	public:

		typedef integer Element;

		inline static bool isUnit (const Element& x) {
			
			return (x == Element(1))  || (x== Element(-1));
		}

		inline static Element& abs(Element& x, const Element& a) {
			x= (a>0)? a: -a;
			return x;
		}

		inline static Element abs(const Element& a){
			return (a>0)? a: -a;
		}

		/** compare two elements, a and b
		  * return 1, if a > b
		  * return 0, if a = b;
		  * return -1. if a < b
		  */
		inline long compare (const Element& a, const Element& b) const {
			
			return (a>b)? 1: ((a<b)? -1 : 0);
		}
		
		/** @memo gcd (g, a, b)
		 *  return g = gcd (a, b)
		 */
		inline static Element& gcd (Element& g, const Element& a, const Element& b) {
				
			Element  u, v, q, r;
			u = a; v = b;

			if (u < 0) {	  
				u = -u;				
			}
 
			if (v < 0) { 	 
				v = -v;
			} 	
 
			while (v != 0) {
				q = u/v;
				r = u -q*v;
				u = v;
				v = r;
			}
 
			g = u;

			return g;
		}
	
		/** @memo gcding (g, b)
		 *  return g = gcd (g, b)
		 */
		inline static Element& gcdin (Element& g, const Element& b) {
			
			gcd(g, g, b);

			return g;
		}

		/** @memo xgcd (g, s, t, a, b)
		 *  g = gcd(a, b) = a*s + b*t.
		 *  The coefficients s and t are defined according to the standard
		 *  Euclidean algorithm applied to |a| and |b|, with the signs then
		 *  adjusted according to the signs of a and b.
		 */
		inline static Element& xgcd (Element& g, Element& s, Element& t, const Element& a, const Element& b){
			Element  u, v, u0, v0, u1, v1, u2, v2, q, r;
 
			int aneg = 0, bneg = 0;
			u = a; v = b;
			if (u < 0) {	  
				u = -u;
				aneg = 1;
			}
 
			if (v < 0) {	 
				v = -v;
				bneg = 1;
			}
 
			u1 = 1; v1 = 0;
			u2 = 0; v2 = 1;

 
			while (v != 0) {
				q = u / v;
				r = u -q*v;
				u = v;
				v = r;
				u0 = u2;
				v0 = v2;
				u2 =  u1 - q*u2;
				v2 = v1- q*v2;
				u1 = u0;
				v1 = v0;
			}
 
			if (aneg)
				u1 = -u1;
 
			if (bneg)
				v1 = -v1;
 
			g = u;
			s = u1;
			t = v1;

			return g;
		}

		/** @memo lcm (c, a, b)
		 *  c = lcm (a, b)
		 */
		inline static Element& lcm (Element& c, const Element& a, const Element& b) {
			
			if ((a==Element(0)) || (b==Element(0))) return c = Element(0);
			
			else {
				Element g;
			
				gcd (g, a, b);
				
				c= a*b;
				c /= g;

				c=abs (c);
			
				return c;
			}
		}
		
		/** @memo lcmin (l, b)
		 *  l = lcm (l, b)
		 */
		inline static Element& lcmin (Element& l, const Element& b) {

			if ((l==Element(0)) || (b==Element(0))) return l = Element(0);
			
			else {
				Element g;
			
				gcd (g, l, b);
				
				l*= b;
				l/= g;

				l=abs (l);
			
				return l;
			}
	
		}


		inline static long reconstructRational (Element& a, Element& b, const Element& x, const Element& m, 
							const Element& a_bound, const Element& b_bound) {
			
			Element  u, v, u0, u1, u2, q, r;
			//cerr<<"approximation: "<<x<<endl<<"base: "<<m<<endl<<"num bound: "<<a_bound<<endl<<"den bound: "<<b_bound<<endl; ;
	      
			u1 = 0; 
			u2 = 1; 
			u = m; v = x;
	      
			while ((v != 0) && ( v > a_bound)) {
				q = u / v;
				r = u -q*v;
				u = v;
				v = r;
				u0 = u2;	 
				u2 =  u1 - q*u2;	 
				u1 = u0;	
			}
	
			if (u2 < Element(0)) { u2= -u2; v=-v;}
			a = v;
			b = u2;

			//cerr<<"rational: "<<a<<"/"<<b<<endl;
	
			return  (b > b_bound)? 0: 1;	

		}


		/** @memo quo (q, x, y)
		 *  q = floor (x/y);
		 */
		inline static Element& quo (Element& q, const Element& a, const Element& b) {
			return  q = a/b;
		}
      
		/** @memo rem (r, a, b)
		 *  r = remindar of  a / b
		 */
		inline static Element& rem (Element& r, const Element& a, const Element& b)  {
			Element q;
			return r= a - quo(q,a,b)*b  ;
		}	

		/** @memo quoin (a, b)
		 *  a = quotient (a, b)
		 */
		inline static Element& quoin (Element& a, const Element& b)  {
			return quo(a,a,b);
		}

		/** @memo quoin (a, b)
		 *  a = quotient (a, b)
		 */
		inline static Element& remin (Element& a, const Element& b)  {
			return rem(a,a,b);
		}

		
		/** @memo quoRem (q, r, a, b)				
		 * q = [a/b], r = a - b*q
		 * |r| < |b|, and if r != 0, sign(r) = sign(b)
		 */
		inline static void quoRem (Element& q, Element& r, const Element& a, const Element& b) {
			quo(q,a,b);
			r = a - q*b;
		}

		/** @memo isDivisor (a, b)
		 *  Test if b | a.
		 */
		inline static bool isDivisor (const Element& a, const Element& b) {
			Element r;
			return rem(r,a,b)==Element(0);
		}

		// some specializations and conversions
		double& convert(double& x, const Element& y) const
		{ return x= (double)y;}

		Element& init(Element& x, const double& y) const 
		{ return x=Element(y);}
      
		integer& convert(integer& x, const Element& y) const
		{ return x=y;}
      
		Element& init(Element& x, const integer& y) const 
		{ return x=y;}


	}; //end of class PID_integer

	template<>
	class ClassifyRing<PID_integer> {
		typedef RingCategories::IntegerTag categoryTag;
	};
	template<>
	std::ostream &UnparametricField<integer>::write (std::ostream &os) const
	{ return os << "unparam<integer>"; }

} //end of namespace LinBox
#endif
