/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/ntl-zz_p.h
 * Copyright (C) 1999-2002 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 */

#ifndef __FIELD_NTL_zz_p_H
#define __FIELD_NTL_zz_p_H

#include <NTL/lzz_p.h>
#include <NTL/ZZ.h>

#include <time.h>
#include "linbox/field/unparametric.h"
#include "linbox/randiter/unparametric.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
  
	/** @name class zz\_p.
	 * Arbitrary precision integers modulus a positive integer.
	 * While NTL allows any integer to serve as the modulus, only prime
	 * moduli yield fields.  Therefore, while arthmetic operations may be
	 * valid for any modulus, only prime moduli are supported in this
	 * implementation.  The primality of the modulus will not be checked, so
	 * it is the programmer's responsibility to supply a prime modulus.
	 * These specializations allow the \Ref{UnparametricField} template class to be
	 * used to wrap NTL's {\tt zz\_p} class as a LinBox field.
	 */
	//@{

	typedef UnparametricField<NTL::zz_p> NTL_zz_p_Field; 

	UnparametricField<NTL::zz_p>::UnparametricField(integer q, size_t e)
	{    
		if(q==0) q=65521;//set default value to 65521
		NTL::zz_p::init(q); // it's an error if q not prime, e not 1.
	}

	/** Initialization of field element from an integer.
	 * Behaves like C++ allocator construct.
	 * This function assumes the output field element x has already been
	 * constructed, but that it is not already initialized.
	 * For now, this is done by converting the integer type to a C++
	 * long and then to the element type through the use of static cast and
	 * NTL's {\tt to\_zz\_p} function.
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
		{ return static_cast<bool>(NTL::IsZero(x)); }

	/** One equality.
	 * Test if field element is equal to one.
	 * This function assumes the field element has already been
	 * constructed and initialized.
	 * In this specialization, NTL's IsOne function is called.
	 * @return boolean true if equals one, false if not.
	 * @param  x field element.
	 */
	template <> bool UnparametricField<NTL::zz_p>::isOne(const NTL::zz_p& x) const
		{ return static_cast<bool>(NTL::IsOne(x)); }

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
	template <> std::ostream& UnparametricField<NTL::zz_p>::write(std::ostream& os) const 
		{ 
			return os << "unparamterized field NTL::zz_p with p = " 
				  << NTL::zz_p::modulus(); 
		}

	/// Constructor for random field element generator
	template <> UnparametricRandIter<NTL::zz_p>::UnparametricRandIter<NTL::zz_p>(
				const UnparametricField<NTL::zz_p>& F, 
				const integer& size, 
				const integer& seed
				)
			: _size(size), _seed(seed)
	{
		if (_seed == integer(0)) _seed = integer(time(NULL));
		
		integer cardinality; 
		F.cardinality(cardinality);
		if (_size > cardinality)
			_size = 0;

#ifdef TRACE
		cout << "created random generator with size " << _size 
   << " and seed " << _seed << endl;
#endif // TRACE
		
		// Seed random number generator
		NTL::SetSeed(NTL::to_ZZ(static_cast<long>(_seed)));
	}

	/// Random field element creator.
	template <> NTL::zz_p& UnparametricRandIter<NTL::zz_p>::random(NTL::zz_p& x)
//		{ return x = static_cast<long>((double(rand())/RAND_MAX)*double(_size)); }
		{
		       if (_size == 0)
		       	       return x = NTL::random_zz_p(); 
		       else
			       return x = NTL::to_zz_p(NTL::RandomBnd(static_cast<long>(_size)));
		}

  

	//@} 
} // namespace LinBox

#endif // __FIELD_NTL_zz_p_H
