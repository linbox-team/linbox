/* -*- mode: c; style: linux -*- */

/* linbox/field/ntl.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __FIELD_NTL_H
#define __FIELD_NTL_H

#include <NTL/RR.h>
#include <NTL/ZZ_p.h>
#include <NTL/lzz_p.h>

#include "linbox/field/unparametric.h"
#include "linbox/randiter/unparametric.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
  
	/** @name NTL Fields.
	 * Partial template instantiantions for wrapping objects from
	 * \URL[Victor Shoup]{http://www.shoup.net/index.html}'s
	 * package \URL[NTL]{http://www.shoup.net/ntl/index.html} 5.0b
	 * objects to comply with the common object interface for
	 * \Ref{LinBox} fields.
	 *
	 * See \URL{http://www.shoup.net/ntl/doc/tour.html} for a tour of the
	 * \URL[NTL]{http://www.shoup.net/ntl/index.html} library.
	 */
	//@{

	/** @name class RR.
	 * Arbitrary precision floating point numbers.
	 * These specializations allow the \Ref{UnparametricField} template class to be
	 * used to wrap NTL's RR class as a LinBox field.
	 */
	//@{

	/** Initialization of field element from an integer.
	 * Behaves like C++ allocator construct.
	 * This function assumes the output field element x has already been
	 * constructed, but that it is not already initialized.
	 * For now, this is done by converting the integer type to a C++
	 * long and then to the element type through the use of static cast and
	 * NTL's to_RR function.
	 * This, of course, assumes such static casts are possible.
	 * This function should be changed in the future to avoid using long.
	 * @return reference to field element.
	 * @param x field element to contain output (reference returned).
	 * @param y integer.
	 */
	template <>
		NTL::RR& UnparametricField<NTL::RR>::init(NTL::RR& x, const integer& y) const
		{ return x = NTL::to_RR(static_cast<const long&>(y)); }

	/** Conversion of field element to an integer.
	 * This function assumes the output field element x has already been
	 * constructed, but that it is not already initialized.
	 * For now, this is done by converting the element type to a C++
	 * long and then to the integer type through the use of static cast and
	 * NTL's to_long function.
	 * This, of course, assumes such static casts are possible.
	 * This function should be changed in the future to avoid using long.
	 * @return reference to integer.
	 * @param x reference to integer to contain output (reference returned).
	 * @param y constant reference to field element.
	 */
	template <>
		integer& UnparametricField<NTL::RR>::convert(integer& x, const NTL::RR& y) const
		{ return x = static_cast<integer>(to_long(y)); }

	/** Multiplicative Inverse.
	 * x = 1 / y
	 * This function assumes both field elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 * @param  y field element.
	 */
	template <> 
		NTL::RR& UnparametricField<NTL::RR>::inv(NTL::RR& x, const NTL::RR& y) const
		{ return x = NTL::inv(y); }
 
	/** Zero equality.
	 * Test if field element is equal to zero.
	 * This function assumes the field element has already been
	 * constructed and initialized.
	 * In this specialization, NTL's IsZero function is called.
	 * @return boolean true if equals zero, false if not.
	 * @param  x field element.
	 */
	template <> bool UnparametricField<NTL::RR>::isZero(const NTL::RR& x) const
		{ return static_cast<bool>(IsZero(x)); }

	/** One equality.
	 * Test if field element is equal to one.
	 * This function assumes the field element has already been
	 * constructed and initialized.
	 * In this specialization, NTL's IsOne function is called.
	 * @return boolean true if equals one, false if not.
	 * @param  x field element.
	 */
	template <> bool UnparametricField<NTL::RR>::isOne(const NTL::RR& x) const
		{ return static_cast<bool>(IsOne(x)); }

	/** Inplace Multiplicative Inverse.
	 * x = 1 / x
	 * This function assumes both field elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 */
	template <> NTL::RR& UnparametricField<NTL::RR>::invin(NTL::RR& x) const
		{ return x = NTL::inv(x); }

	/** Print field.
	 * @return output stream to which field is written.
	 * @param  os  output stream to which field is written.
	 */
	template <> ostream& UnparametricField<NTL::RR>::write(ostream& os) const 
		{ return os << "unparamterized field NTL::RR"; }

	/** Random field element creator.
	 * This returns a random field element from the information supplied
	 * at the creation of the generator.
	 * This generator uses the built-in C++ random number generator instead of
	 * NTL's random function because the NTL function does not allow as much
	 * control over the sampling size as the generic LinBox template.  This 
	 * specialization is included only to allow conversion to an NTL 
	 * object.
	 * @return random field element
	 */
	template <> NTL::RR& UnparametricRandIter<NTL::RR>::random(NTL::RR &elt)
		{
			double temp;
    
			// Create new random elements
			if (_size == 0)
				elt = rand();
			else
				elt = static_cast<long>((double(rand())/RAND_MAX)*double(_size));

			elt = temp;

#ifdef TRACE
			cout << "random double = " << temp 
			     << "    random NTL::RR = " << elt << endl;
#endif // TRACE

			return elt;
    
		} // element& operator() (void)

	//@} // class RR

	/** @name class ZZ_p.
	 * Arbitrary precision integers modulus a positive integer.
	 * While NTL allows any integer to serve as the modulus, only prime
	 * moduli yield fields.  Therefore, while arthmetic operations may be
	 * valid for any modulus, only prime moduli are supported in this
	 * implementation.  The primality of the modulus will not be checked, so
	 * it is the programmer's responsibility to supply a prime modulus.
	 * These specializations allow the \Ref{UnparametricField} template class to be
	 * used to wrap NTL's ZZ_p class as a LinBox field.
	 */
	//@{

	/** Initialization of field element from an integer.
	 * Behaves like C++ allocator construct.
	 * This function assumes the output field element x has already been
	 * constructed, but that it is not already initialized.
	 * For now, this is done by converting the integer type to a C++
	 * long and then to the element type through the use of static cast and
	 * NTL's to_ZZ_p function.
	 * This, of course, assumes such static casts are possible.
	 * This function should be changed in the future to avoid using long.
	 * @return reference to field element.
	 * @param x field element to contain output (reference returned).
	 * @param y integer.
	 */
	template <>
		NTL::ZZ_p& UnparametricField<NTL::ZZ_p>::init(NTL::ZZ_p& x, const integer& y) const
		{ return x = NTL::to_ZZ_p(static_cast<const long&>(y)); }

	/** Conversion of field element to an integer.
	 * This function assumes the output field element x has already been
	 * constructed, but that it is not already initialized.
	 * For now, this is done by converting the element type to a C++
	 * long and then to the integer type through the use of static cast and
	 * NTL's to_long function.
	 * This, of course, assumes such static casts are possible.
	 * This function should be changed in the future to avoid using long.
	 * @return reference to integer.
	 * @param x reference to integer to contain output (reference returned).
	 * @param y constant reference to field element.
	 */
	template <>
		integer& UnparametricField<NTL::ZZ_p>::convert(integer& x, const NTL::ZZ_p& y) const
		{ return x = static_cast<integer>(to_long(rep(y))); }

	/** Cardinality.
	 * Return integer representing cardinality of the field.
	 * Returns the modulus of the field, which should be prime.
	 * @return integer representing cardinality of the field
	 */
	template <> 
		integer& UnparametricField<NTL::ZZ_p>::cardinality(integer& c) const
		{ return c = static_cast<integer>(to_long(NTL::ZZ_p::modulus())); }

	/** Characteristic.
	 * Return integer representing characteristic of the field.
	 * Returns the modulus of the field, which should be prime.
	 * @return integer representing characteristic of the field.
	 */
	template <> 
		integer& UnparametricField<NTL::ZZ_p>::characteristic(integer& c) const
		//FIXME we shouldn't go thru long here as p may be larger than that.
		// check if NTL has cast ZZp to gmp integers.
		{ return c = static_cast<integer>(to_long(NTL::ZZ_p::modulus())); }

	/** Multiplicative Inverse.
	 * x = 1 / y
	 * This function assumes both field elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 * @param  y field element.
	 */
	template <> NTL::ZZ_p& 
		UnparametricField<NTL::ZZ_p>::inv(NTL::ZZ_p& x, const NTL::ZZ_p& y) const
		{ return x = NTL::inv(y); }
 
	/** Zero equality.
	 * Test if field element is equal to zero.
	 * This function assumes the field element has already been
	 * constructed and initialized.
	 * In this specialization, NTL's IsZero function is called.
	 * @return boolean true if equals zero, false if not.
	 * @param  x field element.
	 */
	template <> bool UnparametricField<NTL::ZZ_p>::isZero(const NTL::ZZ_p& x) const
		{ return static_cast<bool>(IsZero(x)); }

	/** One equality.
	 * Test if field element is equal to one.
	 * This function assumes the field element has already been
	 * constructed and initialized.
	 * In this specialization, NTL's IsOne function is called.
	 * @return boolean true if equals one, false if not.
	 * @param  x field element.
	 */
	template <> bool UnparametricField<NTL::ZZ_p>::isOne(const NTL::ZZ_p& x) const
		{ return static_cast<bool>(IsOne(x)); }

	/** Inplace Multiplicative Inverse.
	 * x = 1 / x
	 * This function assumes both field elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 */
	template <> NTL::ZZ_p& UnparametricField<NTL::ZZ_p>::invin(NTL::ZZ_p& x) const
		{ return x = NTL::inv(x); }

	/** Print field.
	 * @return output stream to which field is written.
	 * @param  os  output stream to which field is written.
	 */
	template <> ostream& UnparametricField<NTL::ZZ_p>::write(ostream& os) const 
		{ 
			return os << "unparamterized field NTL::ZZ_p with p = " 
				  << NTL::ZZ_p::modulus(); 
		}

	/** Random field element creator.
	 * This returns a random field element from the information supplied
	 * at the creation of the generator.
	 * This generator uses the built-in C++ random number generator instead of
	 * NTL's random function because the NTL function does not allow as much
	 * control over the sampling size as the generic LinBox template.  This 
	 * specialization is included only to allow conversion to an NTL 
	 * object.
	 * @return random field element
	 */
	template <> NTL::ZZ_p& UnparametricRandIter<NTL::ZZ_p>::random(NTL::ZZ_p &elt)
		{
			long temp;
    
			// Create new random elements
			temp = static_cast<long>((double(rand())/RAND_MAX)*double(_size));

			elt = temp;

#ifdef TRACE
			cout << "random long = " << temp 
			     << "    random NTL::ZZ_p = " << elt << endl;
#endif // TRACE

			return elt;
    
		} // element& operator() (void)

	//@} // class ZZ_p

	/** @name class zz_p.
	 * Arbitrary precision integers modulus a positive integer.
	 * While NTL allows any integer to serve as the modulus, only prime
	 * moduli yield fields.  Therefore, while arthmetic operations may be
	 * valid for any modulus, only prime moduli are supported in this
	 * implementation.  The primality of the modulus will not be checked, so
	 * it is the programmer's responsibility to supply a prime modulus.
	 * These specializations allow the \Ref{UnparametricField} template class to be
	 * used to wrap NTL's zz_p class as a LinBox field.
	 */
	//@{

	/** Initialization of field element from an integer.
	 * Behaves like C++ allocator construct.
	 * This function assumes the output field element x has already been
	 * constructed, but that it is not already initialized.
	 * For now, this is done by converting the integer type to a C++
	 * long and then to the element type through the use of static cast and
	 * NTL's to_zz_p function.
	 * This, of course, assumes such static casts are possible.
	 * This function should be changed in the future to avoid using long.
	 * @return reference to field element.
	 * @param x field element to contain output (reference returned).
	 * @param y integer.
	 */
	template <>
		NTL::zz_p& UnparametricField<NTL::zz_p>::init(NTL::zz_p& x, const integer& y) const
		{ return x = NTL::to_zz_p(static_cast<const long&>(y)); }

	/** Conversion of field element to an integer.
	 * This function assumes the output field element x has already been
	 * constructed, but that it is not already initialized.
	 * For now, this is done by converting the element type to a C++
	 * long and then to the integer type through the use of static cast and
	 * NTL's to_long function.
	 * This, of course, assumes such static casts are possible.
	 * This function should be changed in the future to avoid using long.
	 * @return reference to integer.
	 * @param x reference to integer to contain output (reference returned).
	 * @param y constant reference to field element.
	 */
	template <>
		integer& UnparametricField<NTL::zz_p>::convert(integer& x, const NTL::zz_p& y) const
		{ return x = static_cast<integer>(rep(y)); }

	/** Cardinality.
	 * Return integer representing cardinality of the field.
	 * Returns the modulus of the field, which should be prime.
	 * @return integer representing cardinality of the field
	 */
	template <> 
		integer& UnparametricField<NTL::zz_p>::cardinality(integer& c) const
		{ return c = static_cast<integer>(NTL::zz_p::modulus()); }

	/** Characteristic.
	 * Return integer representing characteristic of the field.
	 * Returns the modulus of the field, which should be prime.
	 * @return integer representing characteristic of the field.
	 */
	template <> 
		integer& UnparametricField<NTL::zz_p>::characteristic(integer& c) const
		{ return c = static_cast<integer>(NTL::zz_p::modulus()); }

	/** Multiplicative Inverse.
	 * x = 1 / y
	 * This function assumes both field elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 * @param  y field element.
	 */
	template <> NTL::zz_p& 
		UnparametricField<NTL::zz_p>::inv(NTL::zz_p& x, const NTL::zz_p& y) const
		{ return x = NTL::inv(y); }
 
	/** Zero equality.
	 * Test if field element is equal to zero.
	 * This function assumes the field element has already been
	 * constructed and initialized.
	 * In this specialization, NTL's IsZero function is called.
	 * @return boolean true if equals zero, false if not.
	 * @param  x field element.
	 */
	template <> bool UnparametricField<NTL::zz_p>::isZero(const NTL::zz_p& x) const
		{ return static_cast<bool>(IsZero(x)); }

	/** One equality.
	 * Test if field element is equal to one.
	 * This function assumes the field element has already been
	 * constructed and initialized.
	 * In this specialization, NTL's IsOne function is called.
	 * @return boolean true if equals one, false if not.
	 * @param  x field element.
	 */
	template <> bool UnparametricField<NTL::zz_p>::isOne(const NTL::zz_p& x) const
		{ return static_cast<bool>(IsOne(x)); }

	/** Inplace Multiplicative Inverse.
	 * x = 1 / x
	 * This function assumes both field elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 */
	template <> NTL::zz_p& UnparametricField<NTL::zz_p>::invin(NTL::zz_p& x) const
		{ return x = NTL::inv(x); }

	/** Print field.
	 * @return output stream to which field is written.
	 * @param  os  output stream to which field is written.
	 */
	template <> ostream& UnparametricField<NTL::zz_p>::write(ostream& os) const 
		{ 
			return os << "unparamterized field NTL::zz_p with p = " 
				  << NTL::zz_p::modulus(); 
		}

	/** Random field element creator.
	 * This returns a random field element from the information supplied
	 * at the creation of the generator.
	 * This generator uses the built-in C++ random number generator instead of
	 * NTL's random function because the NTL function does not allow as much
	 * control over the sampling size as the generic LinBox template.  This 
	 * specialization is included only to allow conversion to an NTL 
	 * object.
	 * @return random field element
	 */
	template <> NTL::zz_p& UnparametricRandIter<NTL::zz_p>::random(NTL::zz_p& x)
		{ return x = static_cast<long>((double(rand())/RAND_MAX)*double(_size)); }

  

	//@} // class zz_p

	//@} NTL Fields

} // namespace LinBox

#endif // __FIELD_NTL_H
