/** -*- mode:C++ -*- */

/** File: ntl-ZZ.h
 *  Author: Zhendong Wan
 */

#ifndef __LINBOX_NTL_ZZ_H__
#define __LINBOX_NTL_ZZ_H__

#include <NTL/ZZ.h>
#include <linbox/integer.h>
#include <iostream>
#include <linbox/util/debug.h>
#include <linbox/randiter/ntl-ZZ.h>

namespace LinBox {
	
	template<class Field>
		class FieldAXPY;
	
	class NTL_ZZ {
		
	public:
		typedef NTL_ZZRandIter RandIter;

		typedef NTL::ZZ Element;

		inline static integer& cardinality (integer& c)  {
			return c = -1;
		}
		
		inline static integer& characteristic (integer& c)  {
			return c = 0;
		}

		static std::ostream& write (std::ostream& out)  {
			return out << "NTL ZZ Ring";
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
		 *   Init from a NTL::ZZ
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
		 *  If z divides y, return x = y / z,
		 *  otherwise, throw an exception
		 */
		inline static Element& div (Element& x, const Element& y, const Element& z) {

			Element q, r;

			NTL::DivRem (q, r, y, z);

			if (NTL::IsZero (r))
				return x = q;

			else
				throw PreconditionFailed(__FUNCTION__,__LINE__,"Div: not dividable");
		}
				
		/** @memo inv (x, y)
		 *  If y is a unit, return x = 1 / y,
		 *  otherwsie, throw an exception
		 */
		inline static Element& inv (Element& x, const Element& y) {

			if ( NTL::IsOne (y)) return x = y;

			else if ( NTL::IsOne (-y)) return x = y;
				
			else 
				throw PreconditionFailed(__FUNCTION__,__LINE__,"Inv: Not invertible");
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
			
			if (NTL::IsOne (x)) return x;
			
			else if (NTL::IsOne (-x)) return x;

			else throw PreconditionFailed(__FUNCTION__,__LINE__,"Div: not dividable");
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
			
			return (NTL::IsOne (x) || NTL::IsOne (-x));
		}
		
		/** @memo gcd (g, a, b)
		 *  return g = gcd (a, b)
		 */
		inline static Element& gcd (Element& g, const Element& a, const Element& b) {
			
			NTL::GCD (g, a, b);

			return g;
		}
	
		/** @memo gcding (g, b)
		 *  return g = gcd (g, b)
		 */
		inline static Element& gcdin (Element& g, const Element& b) {
			
			NTL::GCD (g, g, b);

			return g;
		}

		/** @memo xgcd (g, s, t, a, b)
		 *  g = gcd(a, b) = a*s + b*t.
		 *  The coefficients s and t are defined according to the standard
		 *  Euclidean algorithm applied to |a| and |b|, with the signs then
		 *  adjusted according to the signs of a and b.
		 */
		inline static Element& xgcd (Element& g, Element& s, Element& t, const Element& a, const Element& b){
			
			NTL::XGCD (g,s,t,a,b);

			return g;
		}

		/** @memo lcm (c, a, b)
		 *  c = lcm (a, b)
		 */
		inline static Element& lcm (Element& c, const Element& a, const Element& b) {
			

			if (NTL::IsZero (a) || NTL::IsZero (b)) return c = NTL::ZZ::zero();
			
			else {
				Element g;
			
				NTL::GCD (g, a, b);
				
				NTL::mul (c, a, b);

				c /= g;

				NTL::abs (c, c);
			
				return c;
			}
		}
		
		/** @memo lcmin (l, b)
		 *  l = lcm (l, b)
		 */
		inline static Element& lcmin (Element& l, const Element& b) {

			if (NTL::IsZero (l) || NTL::IsZero (b))
				
				return l = NTL::ZZ::zero();

			else {

				Element g;

				NTL::GCD (g, l, b);

				l *= b;

				l /= g;

				NTL::abs (l, l);

				return l;
			}
		}

				       

				

		// some specail function

		/** @memo sqrt (x, y)
		 *  x = floor ( sqrt(y)).
		 */

		inline static Element& sqrt (Element& x, const Element& y)  {
			
			NTL::SqrRoot(x,y);
			
			return x;
		}
		
		/** @memo  reconstructRational (a, b, x, m, a_bound, b_bound)
		 *  Requires 0 <= x < m, m > 2 * a_bound * b_bound,
		 *  a_bound >= 0, b_bound > 0
		 *   This routine either returns 0, leaving a and b unchanged, 
		 *   or returns 1 and sets a and b so that
		 *  (1) a = b x (mod m),
		 *  (2) |a| <= a_bound, 0 < b <= b_bound, and
		 *  (3) gcd(m, b) = gcd(a, b).
		 */
		
		inline static long reconstructRational (Element& a, Element& b, const Element& x, const Element& m, 
							const Element& a_bound, const Element& b_bound) {
			
			return NTL::ReconstructRational(a,b,x,m,a_bound,b_bound);
		}


		/** @memo quo (q, x, y)
		 *  q = floor (x/y);
		 */
		inline static Element& quo (Element& q, const Element& a, const Element& b) {
			
			NTL::div (q, a, b);

			return q;
		}

		/** @memo rem (r, a, b)
		 *  r = remindar of  a / b
		 */
		inline static Element& rem (Element& r, const Element& a, const Element& b)  {
			
			NTL::rem (r, a, b);
			
			return r;
		}	

		/** @memo quoin (a, b)
		 *  a = quotient (a, b)
		 */
		inline static Element& quoin (Element& a, const Element& b)  {
			
			return a /= b;
			
		}

		/** @memo quoin (a, b)
		 *  a = quotient (a, b)
		 */
		inline static Element& remin (Element& x, const Element& y)  {
			return x %= y;
		}

		
		/** @memo quoRem (q, r, a, b)				
		 * q = [a/b], r = a - b*q
		 * |r| < |b|, and if r != 0, sign(r) = sign(b)
		 */
		inline static void quoRem (Element& q, Element& r, const Element& a, const Element& b) {

			NTL::DivRem(q,r,a,b);
		}

		/** @memo isDivisor (a, b)
		 *  Test if a | b.
		 */
		inline static bool isDivisor (const Element& a, const Element& b) {
			
			if ( NTL::IsZero (a) ) return false;
			
			else if (NTL::IsZero (b)) return true;
			
			else {
				Element r;
				
				NTL::rem (r, b, a);

				return NTL::IsZero (r);
			}
		}
			
	};
		

	template<>
		class FieldAXPY<NTL_ZZ>  {
	public:
		typedef NTL_ZZ Field;
		typedef Field::Element Element;

		/** Constructor.
                 * A faxpy object if constructed from a Field and a field element.
                 * Copies of this objects are stored in the faxpy object.
                 * @param F field F in which arithmetic is done
                 */
                FieldAXPY (const Field &F) : _F (F) { _y = 0; }
 
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
		 * allow optimal multiplication, such as integer * int
                 */
		template<class Element1>
                inline void accumulate (const Element &a, const Element1 &x)
		{ 
			_y += a * x; 
		}
 
                /** Retrieve y
                 *
                 * Performs the delayed modding out if necessary
                 */
                inline Element &get (Element &y) { y = _y; return y; }
 
                /** Assign method.
                 * Stores new field element for arithmetic.
                 * @return reference to self
                 * @param y_init constant reference to element a
                 */
                inline FieldAXPY &assign (const Element& y)
                {
                        _y = y;
                        return *this;
                }
		
		inline void reset() {
			_y = 0;
		}
			
            private:
 
                /// Field in which arithmetic is done
                /// Not sure why it must be mutable, but the compiler complains otherwise
                Field _F;
 
                /// Field element for arithmetic
                Element _y;

	};
}

#endif
