/* linbox/field/ntl-lzz_p.h
 * Copyright (C) 1999-2005 W. J. Turner,
 *               2001 Bradford Hovinen
 * Copyright (C) LinBox
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 * 
 * see COPYING for license information
 *
 */

#ifndef __LINBOX_field_ntl_lzz_p_H
#define __LINBOX_field_ntl_lzz_p_H

#include <NTL/lzz_p.h>
#include <NTL/ZZ.h>

#include <time.h>
#include "linbox/field/unparametric.h"
#include "linbox/randiter/unparametric.h"
#include "linbox/linbox-config.h"
#include <linbox/field/field-traits.h>
#include <linbox/integer.h>

// Namespace in which all LinBox library code resides
namespace LinBox
{
  
  	template <class Ring>
	struct ClassifyRing; 

	class NTL_zz_p;
	
	template<> 
	struct ClassifyRing<NTL_zz_p> {
		typedef RingCategories::ModularTag categoryTag;
	};

/*
	integer& FieldTraits<NTL_zz_p>::maxExponent( integer& i )
		{ return i = integer( "4294967295" ); } // 2^32 - 1
*/

	template<>
	UnparametricField<NTL::zz_p>::UnparametricField(integer q, size_t e)
	{    
		if(q==0) q=65521;//set default value to 65521
		NTL::zz_p::init(q); // it's an error if q not prime, e not 1
	}

	/** Initialization of field element from an integer.
	 * This Uses NTL's <code>to_zz_p</code> function.
	 *
	 * @return reference to field element.
	 * @param x field element to contain output (reference returned).
	 * @param y integer.
	 */
        template <>
        NTL::zz_p& UnparametricField<NTL::zz_p>::init(NTL::zz_p& x, const integer& y) const
                { return x = NTL::to_zz_p(y%NTL::zz_p::modulus()); }
        template <>
        NTL::zz_p& UnparametricField<NTL::zz_p>::init(NTL::zz_p& x, const double& y) const
        {   return x = NTL::to_zz_p((long)(y)%NTL::zz_p::modulus()); }

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
			return os << "unparameterized field NTL::zz_p with p = " 
				  << NTL::zz_p::modulus(); 
		}

	/// Constructor for random field element generator
	template <> 
	UnparametricRandIter<NTL::zz_p>::UnparametricRandIter (const UnparametricField<NTL::zz_p>& F, 
							       const integer& size, 
							       const integer& seed)
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
	template <> NTL::zz_p& UnparametricRandIter<NTL::zz_p>::random(NTL::zz_p& x) const
//		{ return x = static_cast<long>((double(rand())/RAND_MAX)*double(_size)); }
		{
		       if (_size == 0)
		       	       return x = NTL::random_zz_p(); 
		       else
			       return x = NTL::to_zz_p(NTL::RandomBnd(static_cast<long>(_size)));
		}



	/** 
	 * \brief long ints modulo a positive integer.
	 * 
	 * While NTL allows any int to serve as the modulus, only prime
	 * moduli yield fields.  The primality of the modulus will not be checked, so
	 * it is the programmer's responsibility to supply a prime modulus if a field is
	 * wanted.
	 * These specializations allow the \ref{UnparametricField} template class to be
	 * used to wrap NTL's <tt>zz_p</tt> class as a LinBox field.
	 * Uses nice trick for mod p via floating point.
	\ingroup field
	 */

	struct NTL_zz_p: public UnparametricField<NTL::zz_p>
	{
		NTL_zz_p(integer p, size_t e = 1) 
		: UnparametricField<NTL::zz_p>(p, e)
		{}

		NTL::zz_p& init(NTL::zz_p& x, const double& y) const
		{
			double z = fmod(y,NTL::zz_p::modulus());
			if (z > 0) z += 0.5;
			else z -= 0.5;
			return x = NTL::to_zz_p(static_cast<long>(z)); //rounds towards 0
		}

		NTL::zz_p &init (NTL::zz_p &x, const integer &y=0) const {
			NTL::ZZ tmp= NTL::to_ZZ(std::string(y).data());
			return x = NTL::to_zz_p(tmp);
		}


            template <class ANY>
            NTL::zz_p& init(NTL::zz_p& x, const ANY& y) const
	        { return x = NTL::to_zz_p((long)(y)); }

            template <class ANY>
            ANY& convert(ANY& x, const NTL::zz_p& y) const
		{ return x = (ANY)(rep(y)); }
            
	    static inline integer getMaxModulus()
		{ return integer( NTL_SP_BOUND ); }

	    NTL::zz_p& pow( NTL::zz_p& res, const NTL::zz_p& x, long exp ) const
	    {
	    	NTL::power( res, x, exp );
		return res;
	    }

	    NTL::zz_p& powin( NTL::zz_p& x, long exp ) const {
	    	return x = NTL::power(x,exp);
	    }
	};





} // namespace LinBox

#endif // __LINBOX_field_ntl_lzz_p_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
