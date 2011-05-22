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


#ifndef __LINBOX_pid_integer_H
#define __LINBOX_pid_integer_H

#include <limits.h>			    
#include <iostream>
#include <linbox/integer.h>
#include <linbox/field/unparametric.h>
#include <linbox/field/field-traits.h>
#include <linbox/field/gmp-rational.h>


namespace LinBox 
{

	template <class Ring>
	class ClassifyRing;

	/// \ingroup ring
	class PID_integer : public UnparametricField<integer> 
	{

	public:

		//PID_integer(){}
		//PID_integer& operator=(PID_integer& K) { return *this; }

		typedef integer Element;

		inline Element& axpyin (integer &r, const integer& a, const integer& x) const {
			return Integer::axpyin(r,a,x);
		}

                inline Element& axmyin (integer &r, const integer& a, const integer& x) const {
                        return Integer::axmyin(r,a,x);
                }

		inline Element& axpy (integer &r, const integer& a, const integer& x, const integer& y) const {
	                return Integer::axpy(r,a,x,y);//r = ax+y
	        }

		inline  bool isUnit (const Element& x) const { 
			
			return (x == Element(1))  || (x== Element(-1));
		}

		inline  Element& abs(Element& x, const Element& a) const {
			x= (a>0)? a: -a;
			return x;
		}

		inline  Element abs(const Element& a) const {
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
		inline  Element& gcd (Element& g, const Element& a, const Element& b) const {
			return ::gcd(g,a,b);
		}
	
		/** @brief gcdin(g, b)
		 *  return g = gcd (g, b)
		 */
		inline  Element& gcdin (Element& g, const Element& b) const {
			gcd(g, g, b);
			return g;
		}

		/** @brief xgcd (g, s, t, a, b)
		 *  g = gcd(a, b) = a*s + b*t.
		 *  The coefficients s and t are defined according to the standard
		 *  Euclidean algorithm applied to |a| and |b|, with the signs then
		 *  adjusted according to the signs of a and b.
		 */
		inline  Element& xgcd (Element& g, Element& s, Element& t, const Element& a, const Element& b) const {
			return ::gcd(g,a,b,s,t);
		}

		/** @brief lcm (c, a, b)
		 *  c = lcm (a, b)
		 */
		inline  Element& lcm (Element& c, const Element& a, const Element& b) const {
			
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
		
		/** @brief lcmin (l, b)
		 *  l = lcm (l, b)
		 */
		inline  Element& lcmin (Element& l, const Element& b) const {

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
            
                inline  void reconstructRational (Element& a, Element& b, const Element& x, const Element& m) const {
                        RationalReconstruction(a,b, x, m, ::sqrt(m), true, true);
                }
            
                inline  void reconstructRational (Element& a, Element& b, const Element& x, const Element& m, const Element& bound) const {
                        RationalReconstruction(a,b, x, m, bound, true, true);
                }
            
                inline  long reconstructRational (Element& a, Element& b, 
						  const Element& x, const Element& m, 
						  const Element& a_bound, const Element& b_bound) const {
			Element bound = x/b_bound;
			// if (bound>a_bound) std::cerr << "a_bound: " << a_bound << ", x/b_bound: " << bound << std::endl;
			
			RationalReconstruction(a,b,x,m, (bound>a_bound?bound:a_bound), true, false);
			return  (b > b_bound)? 0: 1;	
                }
           
 

		/** @brief quo (q, x, y)
		 *  q = floor (x/y);
		 */
		inline  Element& quo (Element& q, const Element& a, const Element& b) const {
			return  q = a/b;
		}
      
		/** @brief rem (r, a, b)
		 *  r = remindar of  a / b
		 */
		inline  Element& rem (Element& r, const Element& a, const Element& b)  const {
			return Integer::mod(r,a,b);
		}	

		/** @brief quoin (a, b)
		 *  a = quotient (a, b)
		 */
		inline  Element& quoin (Element& a, const Element& b)  const {
			return quo(a,a,b);
		}

		/** @brief quoin (a, b)
		 *  a = quotient (a, b)
		 */
		inline  Element& remin (Element& a, const Element& b)  const {
			return rem(a,a,b);
		}

		
		/** @brief quoRem (q, r, a, b)				
		 * q = [a/b], r = a - b*q
		 * |r| < |b|, and if r != 0, sign(r) = sign(b)
		 */
		inline  void quoRem (Element& q, Element& r, const Element& a, const Element& b) const {
			quo(q,a,b);
			r = a - q*b;
		}

		/** @brief isDivisor (a, b)
		 *  Test if b | a.
		 */
		inline  bool isDivisor (const Element& a, const Element& b) const {
			Element r;
			return rem(r,a,b)==Element(0);
		}

		/** @brief sqrt(x,y)
		 *  x=floor(sqrt(y))
		 */
		inline Element& sqrt(Element& x, const Element& y) const {
			return ::sqrt(x,y);
		}

                inline  Element powtwo(Element& z, const Element& x) const {
                        z = 1;
                        if (x < 0) return z;
						if (x < ULONG_MAX) {
                                z<<=(unsigned long int)x;
                                //cout << "z"<< z;
                                return z;
                        } else {
                                Element n,m;
                                quoRem(n,m,x,(Element)(LONG_MAX-1));
                                for (int i=0; i < n; ++i) {
                                        z <<=(long int)(LONG_MAX-1);
                                }
                                z <= (long int)m;
                                return z;
                        }

                        //for (Element i=0; i < x; ++i) {
                        //      z <<= 1;
                        //}
                        //return z; // BB peut pas !
                }

                inline  Element logtwo(Element& z, const Element& x) const {
                        z = x.bitsize()-1;
	                return z;
			/* 
			if (x<1) return z=-1;
                        z = 0;
                        Element cur = x;
                        cur >>=1;//cout << "cur" << cur;
                        while (cur > 0) {
                                //cout << "cur" << cur;
                                ++z;
                                cur >>=1;
                        }
                        //cout << "z" << z;
                        return z;
			*/
                }



		// some specializations and conversions
		inline double& convert(double& x, const Element& y) const
		{ return x= (double)y;}

		inline Element& init(Element& x, const double& y) const 
		{ return x=Element(y);}
      
		inline integer& convert(integer& x, const Element& y) const
		{ return x=y;}
      
		inline Element& init(Element& x, const integer& y = 0) const 
		{ return x=y;}
		/* 
		 * aniau@astronet.pl 06/2009 initialization form GMPRationalElement
		 */ 
		inline Element& init(Element& x, const GMPRationalElement& q) const
	        {
	                GMPRationalField Q;
	                return Q.convert(x,q);
	        }
            
            inline std::ostream &write (std::ostream &os) const
                { return os << "PID_integer extends unparam<integer>"; }
            
            inline std::ostream &write (std::ostream &os, const Integer& I) const
                { return os << I; }

        protected:
                    // Rational number reconstruction: 
                    // num/den \equiv f modulo m, with |num|<k and 0 < |den| \leq f/k
                    // See [von zur Gathen & Gerhard, Modern Computer Algebra, 
                    //      5.10, Cambridge Univ. Press 1999]
                inline void RationalReconstruction( Element& a, Element& b, 
						    const Element& f, const Element& m, 
						    const Element& k, 
						    bool reduce, bool recursive ) const {
			Element x(f);
                        if (x<0) {
                        	if ((-x)>m)
                            		x %= m;
                        	if (x<0)
                            		x += m;
                    	} else {
                        	if (x>m)
                            		x %= m;
                    	}
			
                        if (x == 0) {
				a = 0;
				b = 1;
                        } else {
				bool res = ratrecon(a,b,x,m,k, reduce, recursive);
				if (recursive)
					for( Element newk = k + 1; (!res) && (newk<f) ; ++newk)
						res = ratrecon(a,b,x,m,newk,reduce, true);
                        }
                }

                // Precondition f is suppposed strictly positive and strictly less than m
                inline  bool ratrecon( Element& num, Element& den, 
                                             const Element& f, const Element& m, 
                                             const Element& k, 
                                             bool reduce, bool recursive ) const {
                    
			//std::cerr << "RatRecon : " << f << " " << m << " " << k << std::endl;                   
                        Element  r0, t0, q, u;
                        r0=m;
                        t0=0;
                        num=f;
                        den=1;
                        while(num>=k)
                        {     
                            q = r0;
                            q /= num;   // r0/num                            
                            u = num;
                            num = r0;  	// num <-- r0
                            r0 = u;	// r0 <-- num
			    Integer::axmyin(num,u,q);
                            if (num == 0) return false;
                            
                            u = den;
                            den = t0;  	// num <-- r0
                            t0 = u;	// r0 <-- num
			    Integer::axmyin(den,u,q);
                        } 

                        if (reduce) {
				// [GG, MCA, 1999] Theorem 5.26
                                				
                                // (ii)
                            Element gg;
                            if (gcd(gg,num,den) != 1) {
                                
                                Element ganum, gar2;
                                for( q = 1, ganum = r0-num, gar2 = r0 ; (ganum < k) && (gar2>=k); ++q ) {
					ganum -= num;
					gar2 -= num;
                                }
                                
                                Integer::axmyin(r0,q,num);
                                Integer::axmyin(t0,q,den);
                                
                                if (t0 < 0) {
                                    num = -r0;
                                    den = -t0;
                                } else {
                                    num = r0;
                                    den = t0;
                                }
                                
//                                if (t0 > m/k) {
				if (den > m/k) {
                                    if (!recursive) 
                                        std::cerr 
                                            << "*** Error *** No rational reconstruction of " 
                                            << f 
                                            << " modulo " 
                                            << m 
                                            << " with denominator <= " 
                                            << (m/k)
                                            << std::endl;   
                                }
                                if (gcd(gg,num,den) != 1) {
                                    if (!recursive) 
                                        std::cerr 
                                            << "*** Error *** There exists no rational reconstruction of " 
                                            << f 
                                            << " modulo " 
                                            << m 
                                            << " with |numerator| < " 
                                            << k
                                            << std::endl
                                            << "*** Error *** But " 
                                            << num
                                            << " = "
                                            << den
                                            << " * "
                                            << f
                                            << " modulo " 
                                            << m 
                                            << std::endl;
                                    return false;
                                }
                            }
                        }
                       // (i)
                       if (den < 0) {
	                       Integer::negin(num);
                               Integer::negin(den);
                       }

// std::cerr << "RatRecon End " << num << "/" << den << std::endl;
                        return true;    
                }
            
	}; //end of class PID_integer

	template<>
	struct ClassifyRing<PID_integer> {
		typedef RingCategories::IntegerTag categoryTag;
	};


	/*
	// Specialization for Homomorphism
	template <class _Target>
	class Hom<PID_integer, _Target> {
	public:
		typedef PID_integer Source;
		typedef _Target Target;
		typedef typename Source::Element SrcElt;
		typedef typename Target::Element Elt;
		
		Hom(const Source& S, const Target& T) : _source (S), _target(T){}

		Elt& image(Elt& t, const SrcElt& s) {
			if  (s.bitsize() > 52 )
				_target.init(t,s);
			else
				_target.init(t, (double)s);
			return t;
		}

		SrcElt& preimage(SrcElt& s, const Elt& t) {
			_source.convert(s,t);
			return s;
		}

		const Source& source() { return _source;}

		const Target& target() { return _target;}
		
	protected:
		double tmp;
		Source _source;
		Target _target;
		
	};
	
	*/
} //end of namespace LinBox

#endif //__LINBOX_pid_integer_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
