/* linbox/field/modular.h
 * Copyright(C) LinBox
 * Written by 
 *    Pierrick Vignard 
 *    Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 *
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_poweroftwomodular_H
#define __LINBOX_poweroftwomodular_H

#include <iostream>

#include "linbox/integer.h"
#include "linbox/ring/gcd32.h"
#include "linbox/ring/gcd64.h"
#include "linbox/randiter/unparametric.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 
	
	/** \brief Ring of elements modulo some power of two
	\ingroup ring
	 *
	 * @param element Element type, e.g. long or integer
	 * @param Intermediate Type to use for intermediate computations. This
	 *                     should be a data type that can support integers
	 *                     twice the length of the maximal modulus used
	 */
  template <class Ints> class PowerOfTwoModular 
	{
	public:
		
		/** Element type
		 */
		typedef Ints Element;

		/** Random iterator generator type.
		 * It must meet the common object interface of random element generators
		 * as given in the the archetype RandIterArchetype.
		 */
	        struct RandIter{
		  typedef Ints Element;

		  RandIter ( const PowerOfTwoModular<Ints>& F, 
			     const integer& size = 0, const integer& seed = 0){
		    if (_seed == integer(0)) _seed = integer(time(NULL));
		    srand(static_cast<long>(_seed));
			
		  }

		  RandIter ( const RandIter& R ):_seed(R._seed){  
		    
		  }

		  RandIter ( void ):_seed(0){  
		    
		  }

		Element& random (Element& x) const
		{
		  return x=rand();
		}
		protected:
		  integer _seed;
		};

		/** @name Object Management
		 */
		//@{
 
		/** Default constructor.
		 */
		PowerOfTwoModular (void) {
			_poweroftwo=sizeof(Ints)<<3;
		}

 
		/** Conversion of field base element to a template class T.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.
		 * @return reference to template class T.
		 * @param x template class T to contain output (reference returned).
		 * @param y constant field base element.
		 */
		integer &convert (integer &x, const Element &y) const
		{ return x = integer (y); }

		/** Initialization of field base element from an integer.
		 * Behaves like C++ allocator construct.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.
		 * This is not a specialization of the template function because
		 * such a specialization is not allowed inside the class declaration.
		 * @return reference to field base element.
		 * @param x field base element to contain output (reference returned).
		 * @param y integer.
		 */
		Element &init (Element &x, const Ints &y = 0) const
		{ 
			return x=y;
		} 


		/** Assignment of one field base element to another.
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &assign (Element &x, const Element &y) const { return x = y; }

		/** Cardinality.
		 * Return integer representing cardinality of the domain.
		 * Returns a non-negative integer for all domains with finite
		 * cardinality, and returns -1 to signify a domain of infinite
		 * cardinality.
		 * @return integer representing cardinality of the domain
		 */
		integer &cardinality (integer &c) const
		{ return c = integer (1) <<_poweroftwo; }

		/** Characteristic.
		 * Return integer representing characteristic of the domain.
		 * Returns a positive integer to all domains with finite characteristic,
		 * and returns 0 to signify a domain of infinite characteristic.
		 * @return integer representing characteristic of the domain.
		 */
		integer &characteristic (integer &c) const
		{ return c = integer (1) <<_poweroftwo; }


		/** poweroftwo */
		int &poweroftwo(int &c)
		{ return c=_poweroftwo; }

		//@} Object Management

		/** @name Arithmetic Operations
		 * x <- y op z; x <- op y
		 * These operations require all elements, including x, to be initialized
		 * before the operation is called.  Uninitialized field base elements will
		 * give undefined results.
		 */
		//@{

		/** Equality of two elements.
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return boolean true if equal, false if not.
		 * @param  x field base element
		 * @param  y field base element
		 */
		bool areEqual (const Element &x, const Element &y) const
		{ return x == y; }

		/** Zero equality.
		 * Test if field base element is equal to zero.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return boolean true if equals zero, false if not.
		 * @param  x field base element.
		 */
		bool isZero (const Element &x) const
		{ return x == 0; }
 
		/** One equality.
		 * Test if field base element is equal to one.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return boolean true if equals one, false if not.
		 * @param  x field base element.
		 */
		bool isOne (const Element &x) const
		{ return x == 1; }

		
		bool isUnit (const Element &x) const
		  { return x&1UL; }

		bool isZeroDivisor ( const Element &x ) const
		  { return !(x&1UL); }

		/** Gcd with 2^_poweroftwo *
		 * Valid for Ints up to 32 bits *
		 * Specialization is required for bigger Ints
		 */
		Element &gcd_poweroftwo (Element &x,const Element &y) const
		{ 
		  return x=GCD2E32(y);
		}

		/** Does x divide y */
		bool doesdivide (const Element &x, const Element &y) const
		{ 
			Element tmp1,tmp2;
			gcd_poweroftwo(tmp1,x),gcd_poweroftwo(tmp2,y);
			return tmp1<=tmp2;
		}

		/** Power of two in x 
		 * Input Element x = 2^n*y where y is odd
		 * Output n
		 */
		int poweroftwoinx (const Element &x) const
		{
			Element tmp;
			gcd_poweroftwo(tmp,x);
			//printf("tmp = %d\n",tmp);
			int n=0;
			while(tmp^1)
				{
					tmp>>=1;
					++n;
					//printf("tmp = %d et n = %d\n",tmp,n);
				}
			return n;
					
		}

		Element &bezout(const Element &x, const Element &y, Element &gcd, Element &u, Element &v) const
		{
			Element v1,v2,v3,t1,t2,t3,q;
			u=1;
			v=0;
			gcd=x;
			v1=0;
			v2=1;
			v3=y;
			while(v3!=0)
				{
					q=gcd/v3;
					t1=u-q*v1;
					t2=v-q*v2;
					t3=gcd-q*v3;
					u=v1;
					v=v2;
					gcd=v3;
					v1=t1;
					v2=t2;
					v3=t3;
				}
			return gcd;
		}

		//@} Arithmetic Operations

		/** @name Input/Output Operations */
		//@{

		/** Print field.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 */
		std::ostream &write (std::ostream &os) const 
		{ return os << "integers mod 2^" << _poweroftwo; }

		/** Read field.
		 * @return input stream from which field is read.
		 * @param  is  input stream from which field is read.
		 */
		std::istream &read (std::istream &is) { return is >> _poweroftwo; }

		/** Print field base element.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return output stream to which field base element is written.
		 * @param  os  output stream to which field base element is written.
		 * @param  x   field base element.
		 */
		std::ostream &write (std::ostream &os, const Element &x) const
		{ return os << x; }
 
		/** Read field base element.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return input stream from which field base element is read.
		 * @param  is  input stream from which field base element is read.
		 * @param  x   field base element.
		 */
		std::istream &read (std::istream &is, Element &x) const
		{
			return is >> x;
		}



		//@}  

		/** @name Arithmetic Operations
		 * x <- y op z; x <- op y
		 * These operations require all elements, including x, to be initialized
		 * before the operation is called.  Uninitialized field base elements will
		 * give undefined results.
		 */
		//@{

		/** Addition.
		 * x = y + z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element &add (Element &x, const Element &y, const Element &z) const
		{
			return x = y + z;
		}
 
		/** Subtraction.
		 * x = y - z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element &sub (Element &x, const Element &y, const Element &z) const
		{ 
			return x = y - z;
		}
 
		/** Multiplication.
		 * x = y * z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element &mul (Element &x, const Element &y, const Element &z) const
		{ return x = (y * z); }
 
		/** Division.
		 * x = y / z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * This fonction assumes that x divides y. That can be verified by using doesdivide(x,y)
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element &div (Element &x, const Element &y, const Element &z) const
		{ 
			int n=poweroftwoinx(z);
			Element tmp;
			inv(tmp, z>>n );
			return mul(x, y>>n , tmp);
		}

		/** Additive Inverse (Negation).
		 * x = - y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &neg (Element &x, const Element &y) const
		{ return x = (~y)+1; }
		//{ return x = 0-y; }
 
		/** Multiplicative Inverse.
		 * x = 1 / y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * This function assumes that y is odd (ie 1/y exists)
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &inv (Element &x, const Element &y) const
		{
			Element u,v;
			neg(v,y>>1);
			u=(y+(v<<2));
			for(unsigned int puiss=2;puiss<_poweroftwo;puiss<<=1)
				{
					v*=v;
					u*=(u*y+(v<<(puiss+1)));
				}
			return x=u;
		}


		/** Multiplicative Inverse 2.
		 * x = 1 / y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &inv2 (Element &x, const Element &y) const
		{
			Element gcd,q,r,u,v,twoexpnmu=1;
			twoexpnmu=twoexpnmu<<(_poweroftwo-1);
			q=twoexpnmu/y;
			r=twoexpnmu-y*q;
			bezout(y,r<<1,gcd,u,v);
			return x=u-((q*v)<<1);
		}

		/** Natural AXPY.
		 * r  = a * x + y
		 * This function assumes all field elements have already been 
		 * constructed and initialized.
		 * @return reference to r.
		 * @param  r field element (reference returned).
		 * @param  a field element.
		 * @param  x field element.
		 * @param  y field element.
		 */
		Element &axpy (Element &r, 
			       const Element &a, 
			       const Element &x, 
			       const Element &y) const
		{ 
			return r = (a * x + y);
		}

		//@} Arithmetic Operations
 
		/** @name Inplace Arithmetic Operations
		 * x <- x op y; x <- op x
		 */
		//@{

		/** Inplace Addition.
		 * x += y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &addin (Element &x, const Element &y) const
		{ 
			return x += y;
		}
 
		/** Inplace Subtraction.
		 * x -= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &subin (Element &x, const Element &y) const
		{
			return x -= y;
		}
 
		/** Inplace Multiplication.
		 * x *= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &mulin (Element &x, const Element &y) const
		{
			return x*=y;
		}
 
		/** Inplace Division.
		 * x /= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &divin (Element &x, const Element &y) const
		{
			Element temp;
			inv (temp, y);
			return mulin (x, temp);
		}
 
		/** Inplace Additive Inverse (Inplace Negation).
		 * x = - x
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 */
		Element &negin (Element &x) const
		{
			return x = (~x)+1;
		}
 
		/** Inplace Multiplicative Inverse.
		 * x = 1 / x
		 * This function assumes the field base elementhas already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 */
		Element &invin (Element &x) const
		{ return inv (x, x); }

		/** Inplace AXPY.
		 * r  += a * x
		 * This function assumes all field elements have already been 
		 * constructed and initialized.
		 * Purely virtual
		 * @return reference to r.
		 * @param  r field element (reference returned).
		 * @param  a field element.
		 * @param  x field element.
		 */
		Element &axpyin (Element &r, const Element &a, const Element &x) const
		{ 
			return r = (r + a * x);
		}

		//@} Inplace Arithmetic Operations


		//@}

	protected:

		/// Private (non-static) element for modulus
		Element _poweroftwo;

	}; // class PowerOfTwoModular
	
#if 0
	/* Specialization of gcd_poweroftwo for Int64
	*/
	template<> 
	PowerOfTwoModular<int64>::Element&  
	PowerOfTwoModular<int64>::gcd_poweroftwo (Element &x,const Element &y) const 
	{  
		return x=GCD2E64(y); 
	}  
#endif

} // namespace LinBox

// #include "linbox/field/modular.inl"
// #include "linbox/randiter/modular.h"

#endif // __LINBOX_poweroftwomodular_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
