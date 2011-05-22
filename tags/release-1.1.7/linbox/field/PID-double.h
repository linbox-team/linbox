/* linbox/field/PID-double.h
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




#ifndef __LINBOX_pid_double_H
#define __LINBOX_pid_double_H

#include <linbox/field/unparametric.h>
#include <linbox/field/field-traits.h>

namespace LinBox 
{


	/// \ingroup ring
	class PID_double : public UnparametricField<double> 
	{

	public:

		typedef double Element;

		inline static bool isUnit (const Element& x) {
			
			return (x == double(1))  || (x== double(-1));
		}

		inline static Element& abs (Element& x, const Element& a) {
			x= (a>0)? a: -a;
			return x;
		}

		inline static Element abs (const Element& a){
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
		

		/** @brief gcd (g, a, b)
		 *  return g = gcd (a, b)
		 */
		inline static Element& gcd (Element& g, const Element& a, const Element& b) {
				
			double  u, v, q, r;
			u = a; v = b;

			if (u < 0) {	  
				u = -u;				
			}
 
			if (v < 0) { 	 
				v = -v;
			} 	
 
			while (v != 0) {
				q = floor(u/v);
				r = u -q*v;
				u = v;
				v = r;
			}
 
			g = u;

			return g;
		}
	
		/** @brief gcding (g, b)
		 *  return g = gcd (g, b)
		 */
		inline static Element& gcdin (Element& g, const Element& b) {
			
			gcd(g, g, b);

			return g;
		}

		/** @brief xgcd (g, s, t, a, b)
		 *  g = gcd(a, b) = a*s + b*t.
		 *  The coefficients s and t are defined according to the standard
		 *  Euclidean algorithm applied to |a| and |b|, with the signs then
		 *  adjusted according to the signs of a and b.
		 */
		inline static Element& xgcd (Element& g, Element& s, Element& t, const Element& a, const Element& b){
			double  u, v, u0, v0, u1, v1, u2, v2, q, r;
 
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
				q = floor(u / v);
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

		/** @brief lcm (c, a, b)
		 *  c = lcm (a, b)
		 */
		inline static Element& lcm (Element& c, const Element& a, const Element& b) {
			
			if ((a==0.) || (b==0.)) return c = 0.;
			
			else {
				Element g;
			
				gcd (g, a, b);
				
				c= a*b;
				c /= g;

				if (c<0) c = -c; //no abs for double, dpritcha
				// c=abs (c);
			
				return c;
			}
		}
		
		/** @brief lcmin (l, b)
		 *  l = lcm (l, b)
		 */
		inline static Element& lcmin (Element& l, const Element& b) {

			if ((l==0.) || (b==0.)) return l = 0.;
			
			else {
				Element g;
			
				gcd (g, l, b);
				
				l*= b;
				l/= g;
	  
				//l=abs (l);
				if (l < 0) l = -l;
				return l;
			}
	
		}


		inline static long reconstructRational (Element& a, Element& b, const Element& x, const Element& m, 
							const Element& a_bound, const Element& b_bound) {
			
			double  u, v, u0, u1, u2, q, r;
	
			u1 = 0; 
			u2 = 1; 
			u = m; v = x;
 
			while ((v != 0) && ( v > a_bound)) {
				q = floor(u / v);
				r = u -q*v;
				u = v;
				v = r;
				u0 = u2;	 
				u2 =  u1 - q*u2;	 
				u1 = u0;	
			}
	
			if (u2 < 0.) { u2= -u2; v=-v;}
			a = v;
			b = u2;
	
			return  (b > b_bound)? 0: 1;	

		}


		/** @brief quo (q, x, y)
		 *  q = floor (x/y);
		 */
		inline static Element& quo (Element& q, const Element& a, const Element& b) {
			return  q = floor (a/b);
		}
      
		/** @brief rem (r, a, b)
		 *  r = remindar of  a / b
		 */
		inline static Element& rem (Element& r, const Element& a, const Element& b)  {
			Element q;
			return r= a - quo(q,a,b)*b  ;
		}	

		/** @brief quoin (a, b)
		 *  a = quotient (a, b)
		 */
		inline static Element& quoin (Element& a, const Element& b)  {
			return quo(a,a,b);
		}

		/** @brief quoin (a, b)
		 *  a = quotient (a, b)
		 */
		inline static Element& remin (Element& a, const Element& b)  {
			return rem(a,a,b);
		}

		
		/** @brief quoRem (q, r, a, b)				
		 * q = [a/b], r = a - b*q
		 * |r| < |b|, and if r != 0, sign(r) = sign(b)
		 */
		inline static void quoRem (Element& q, Element& r, const Element& a, const Element& b) {
			quo(q,a,b);
			r = a - q*b;
		}

		/** @brief isDivisor (a, b)
		 *  Test if a | b.
		 */
		inline static bool isDivisor (const Element& a, const Element& b) {
			double r;
			return rem(r,a,b)==0.;
		}

		// some specializations and conversions
		double& convert(double& x, const Element& y) const
		{ return x=y;}

		Element& init(Element& x, const double& y) const 
		{ return x=y;}
      
		integer& convert(integer& x, const Element& y) const
		{ return x=(integer)y;}
      
		Element& init(Element& x, const integer& y) const 
		{ return x=(double)y;}


	}; //end of class PID_double

	template<>
	std::ostream &UnparametricField<double>::write (std::ostream &os) const
	{ return os << "unparam<double>"; }

} //end of namespace LinBox
#endif //__LINBOX_pid_double_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
