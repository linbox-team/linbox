/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/PIR-ntl-ZZ_p.h
 * written by Zhendong Wan
 *
 */

#ifndef __PIR_NTL_ZZ_p_H__
#define __PIR_NTL_ZZ_p_H__

#include <linbox/field/unparametric.h>
#include "linbox-config.h"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

// Namespace in which all LinBox library code resides
namespace LinBox
{
	
	template<class Field>
	class FieldAXPY;
	
	/** @memo extend Wrapper of ZZ_p from NTL.  Add PIR functions
	 */

	class PIR_ntl_ZZ_p : public UnparametricField<NTL::ZZ_p> {
		
	public:
		typedef NTL::ZZ_p Element;

		PIR_ntl_ZZ_p (const NTL::ZZ& d) {

			NTL::ZZ_p::init(d);

		}
		

		inline static integer& cardinality (integer& c);
			
		inline static NTL::ZZ& cardinality (NTL::ZZ& c) {

			return c = NTL::ZZ_p::modulus();
		}
		
		inline static integer& characteristic (integer& c);


		static std::ostream& write (std::ostream& out)  {
			return out << "PIR_NTL_ZZ_p Ring";
		}

		std::istream& read (std::istream& in)  {
			return in;
		}
		
		/** @memo Init (x, y)
		 *  Init x from y.
		 */
		template<class Element2>
			inline static Element& init (Element& x,  const Element2& y) {
			
			NTL::conv (x, y);
			
			return x;
		}

		/** @memo Specialization of init.
		 *   Init from a NTL::ZZ_p
                 */
                inline static Element& init (Element& x, const Element& y) {
	
			x = y;

			return x;
		}

		/** @memo Specialization of init.
		 *  I don't  know how to init from integer.
		 */
		inline static Element& init (Element& x, const integer& y) ;
		
		/** @memo Convert (x, y).
		 *  Convert y to an Element.
		 */
		static integer& convert (integer& x, const Element& y)  ;
		
		/** @memo Assign (x, y);
		 *  x = y.
		 */
		inline static Element&  assign (Element& x, const Element& y)  {
			return x = y;
		}

		/** @memo areEqual (x, y)
		 *  Test if x == y
		 */
		inline static bool areEqual (const Element& x ,const Element& y)  {
			return x == y;
		}

		/** @memo isZero (x)
		 *  Test if x == 0
		 */
		inline static bool isZero (const Element& x)  {
			
			return NTL::IsZero (x);
		}

		/** @memo isOne (x)
		 *  Test if x == 1
		 */
		inline static bool isOne (const Element& x)  {

			return NTL::IsOne (x);
		}
								
		// arithmetic
		
		/** @memo add (x, y, z)
		 *  return x = y + z
		 */
		inline static Element& add (Element& x, const Element& y, const Element& z)  {

			NTL::add (x, y, z);

			return x;
		}

		/** @memo sub (x, y, z)
		 *  return x = y - z
		 */
		inline static Element& sub (Element& x, const Element& y, const Element& z)  {			
			
			NTL::sub (x, y, z);

			return x;
		}
			      
		/** @memo mul (x, y, z)
		 *  return x = y * z
		 */
		template <class Int>
			inline static Element& mul (Element& x, const Element& y, const Int& z)  {
			
			NTL::mul (x, y, z);

			return x;
		}

		/** @memo div (x, y, z)
		 *  If exists a, such that a * z =y,
		 *  return x = one of them.
		 *  Otherwise, throw an exception
		 */
		inline static Element& div (Element& x, const Element& y, const Element& z) {

			NTL::ZZ g, s, t;

			NTL::XGCD (g, s, t, NTL::rep(z), NTL::ZZ_p::modulus());
			
			NTL::ZZ q, r;

			NTL::DivRem (q, r, NTL::rep(y), g);

			if (NTL::IsZero (r)) {

				Element tmp1, tmp2;

				NTL::conv (tmp1, s);

				NTL::conv (tmp2, q);

				NTL::mul (x, tmp1, tmp2);
				
			} 
					
			else
				throw PreconditionFailed(__FUNCTION__,__LINE__,"Div: not dividable");


			return x;
				
		}
				
		/** @memo inv (x, y)
		 *  If y is a unit, return x = 1 / y,
		 *  otherwsie, throw an exception
		 */
		inline static Element& inv (Element& x, const Element& y) {

			NTL::inv (x, y);

			return x;
		}

		/** @memo neg (x, y)
		 *  return x = -y;
		 */
		inline static Element& neg (Element& x, const Element& y)  {
			
			NTL::negate (x, y);

			return x;
		}


		/** @memo axpy (r, a, x, y)
		 *  return r = a x + y
		 */

		template <class Int>
			inline static Element& axpy (Element& r, const Element& a, const Int& x, const Element& y)  {

			NTL::mul (r, a, x);

			return r += y;
		}


		// inplace operator
		
		/** @memo addin (x, y)
		 *  return x += y;
		 */
		inline static Element& addin (Element& x, const Element& y) {
			
			return x += y;
		}
		
		/** @memo subin (x, y)
		 *  return x -= y;
		 */
		inline static Element& subin (Element& x, const Element& y)  {
			
			return x -= y;
		}

		/** @memo mulin (x, y)
		 *  return x *= y;
		 */
		template<class Int>
			inline static Element& mulin (Element& x, const Int& y)  {
			
			return x *= y;
		}

		/** @memo divin (x, y)
		 *  If y divides x, return x /= y,
		 *  otherwise throw an exception
		 */
		inline static Element& divin (Element& x, const Element& y) {
			
			div (x, x, y);

			return x;
		}

		/** @memo invin (x)
		 *  If x is a unit, x = 1 / x,
		 *  otherwise, throw an exception.
		 */
		inline static Element& invin (Element& x) {
			
			return x = NTL::inv(x);
		}				
		
		/** @memo negin (x)
		 *  return x = -x;
		 */
		inline static Element& negin (Element& x)  {			

			NTL::negate (x, x);

			return x;
		}

		/** @memo axpyin (r, a, x)
		 *  return r += a x
		 */
		template <class Int>
			inline static Element& axpyin (Element& r, const Element& a, const Int& x)  {

			return r += a * x;
		}

	
		// IO

		/** @memo write (out, y)
		 *  out << y;
		 */
		static std::ostream& write(std::ostream& out,const Element& y)  {

			out << y;
			
			return out;
		}


		/** @memo read (in, x)
		 *  read x from istream in
		 */
		static std::istream& read(std::istream& in, Element& x) {
			
			return in >> x;
		}


		/** some PIR function
		 */

		/** @memo isUnit (x)
		 *  Test if x is a unit.
		 */
		inline static bool isUnit (const Element& x) {
			
			NTL::ZZ g;

			NTL::GCD (g, NTL::rep(x), NTL::ZZ_p::modulus());

			return NTL::IsOne (g);
		}
		
		/** @memo gcd (g, a, b)
		 *  return g = gcd (a, b)
		 */
		inline static Element& gcd (Element& g, const Element& a, const Element& b) {
			
			NTL::ZZ d;
			
			NTL::GCD (d, NTL::rep(a), NTL::rep(b));

			NTL::conv (g, d);

			return g;
		}
	
		/** @memo gcding (g, b)
		 *  return g = gcd (g, b)
		 */
		inline static Element& gcdin (Element& g, const Element& b) {
			
			gcd (g, g, b);

			return g;
		}

		/** @memo xgcd (g, s, t, a, b)
		 *  g = gcd(a, b) = a*s + b*t.
		 *  and gcd (s, t) is a unit.
		 */
		inline static Element& xgcd (Element& g, Element& s, Element& t, const Element& a, const Element& b){

			NTL::ZZ g1, s1, t1;
			
			NTL::XGCD (g1, s1, t1, NTL::rep(a), NTL::rep(b));

			NTL::conv (g, g1);

			NTL::conv (s, s1);

			NTL::conv (t, t1);

			return g;
		}

		/** @memo dxgcd (g, s, t, a1, b1, a, b)
		 *  g = gcd(a, b) = a*s + b*t.
		 *  and gcd (s, t) is a unit.
		 *  s * a1 + t * b1 = a unit.
		 */
		inline static Element& dxgcd (Element& g, Element& s, Element& t, Element& a1, Element& b1, 
					      const Element& a, const Element& b){

			NTL::ZZ g1, s1, t1, a2, b2;
			
			NTL::XGCD (g1, s1, t1, NTL::rep(a), NTL::rep(b));

			NTL::conv (g, g1);

			NTL::conv (s, s1);

			NTL::conv (t, t1);

			NTL::div (a2, NTL::rep(a), g1);

			NTL::div (b2, NTL::rep(b), g1);

			NTL::conv (a1, a2);

			NTL::conv (b1, b2);

			return g;
		}

		/** @memo isDivisor (a, b)
		 *  Test if a | b.
		 */
		inline static bool isDivisor (const Element& a, const Element& b) {
			
			if ( NTL::IsZero (a) ) return false;
			
			else if (NTL::IsZero (b)) return true;
			
			else {
				NTL::ZZ g, r;

				NTL::GCD (g, NTL::rep(a), NTL::ZZ_p::modulus());
								
				NTL::rem (r, NTL::rep(b), g);

				return NTL::IsZero (r);
			}
		}

		/** @memo normal (a, b);
		 *  a = normalization of b.
		 */
			
		inline static Element& normal (Element& a,  const Element& b) {

			NTL::ZZ a1;

			NTL::GCD (a1, NTL::rep(b), NTL::ZZ_p::modulus());

			NTL::conv (a, a1);
			
			return a;
		}

		/** @memo normalIn (a)
		 */

		inline static Element& normalIn (Element& a) {

		
			NTL::ZZ a1;

			NTL::GCD (a1, NTL::rep(a), NTL::ZZ_p::modulus());

			NTL::conv (a, a1);
			
			return a;	
		}
		    
	};

	
		
	template <>
	class FieldAXPY<PIR_ntl_ZZ_p>  {
	public:
		typedef PIR_ntl_ZZ_p Field;
		typedef Field::Element Element;

		/** Constructor.
                 * A faxpy object if constructed from a Field.
                 * Copies of this objects are stored in the faxpy object.
                 * @param F field F in which arithmetic is done
                 */
                FieldAXPY (const Field &F) : _F (F) { _y = NTL::ZZ::zero(); }
 
                /** Copy constructor.
                 * @param faxpy
                 */
                FieldAXPY (const FieldAXPY<Field> &faxpy) : _F (faxpy._F), _y (faxpy._y) {}
 
                /** Assignment operator
                 * @param faxpy
                 */
                FieldAXPY<Field> &operator = (const FieldAXPY &faxpy)
                        { _y = faxpy._y; return *this; }
 
                /** Add a*x to y
                 * y += a*x.
                 * @param a constant reference to element a
                 * @param x constant reference to element x
		 */
		inline void accumulate (const Element &a, const Element &x)
		{ 
			_y += NTL::rep(a) * NTL::rep(x); 
		}
 
                /** Retrieve y
                 *
                 * Performs the delayed modding out if necessary
                 */
                inline Element &get (Element &y) { NTL::conv (y,  _y); return y; }
 
                /** Assign method.
                 * Stores new field element for arithmetic.
                 * @return reference to self
                 * @param y_init constant reference to element a
                 */
                inline FieldAXPY &assign (const Element& y)
                {
                        _y = NTL::rep(y);
			
                        return *this;
                }
		
		inline void reset() {
			_y = NTL::ZZ::zero();
		}
			
            private:
 
                /// Field in which arithmetic is done
                /// Not sure why it must be mutable, but the compiler complains otherwise
                Field _F;
 
                /// Field element for arithmetic
                NTL::ZZ _y;

	};
}

#endif
