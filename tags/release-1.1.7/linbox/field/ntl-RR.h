/* linbox/field/ntl-RR.h
 * Copyright (C) 1999-2005 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_field_ntl_rr_H
#define __LINBOX_field_ntl_rr_H
#include <NTL/tools.h>

#include <NTL/RR.h>

#include "linbox/field/unparametric.h"
#include "linbox/linbox-config.h"
#include <linbox/field/field-traits.h>

// Namespace in which all LinBox library code resides
namespace LinBox
{
  
	template <class Ring>
	struct ClassifyRing;

	template<class Element>
	struct ClassifyRing<UnparametricField<Element> >;

	template<>
	struct ClassifyRing<UnparametricField<NTL::RR> >{
		typedef RingCategories::ModularTag categoryTag;
	};

	//class NTL_RR: public UnparametricField<NTL::RR>, public FieldInterface{};
	//typedef UnparametricField<NTL::RR> NTL_RR;

	/** @name class RR.
	 \brief
	 * Rational number field.  
         * This field is provided as a convenience in a few places.  
         * Use with caution because expression swell.
         *
	 * This specialization allows the \ref{UnparametricField} template class to be
	 * used to wrap NTL's RR class as a LinBox field.
	\ingroup field
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
	template <>
		NTL::RR& UnparametricField<NTL::RR>::init(NTL::RR& x, const double& y) const
		{ return x = NTL::to_RR((long)(y)); }

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
	template <> std::ostream& UnparametricField<NTL::RR>::write(std::ostream& os) const 
		{ return os << "unparameterized field NTL::RR"; }


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
	template <> NTL::RR& UnparametricRandIter<NTL::RR>::random(NTL::RR &elt) const
		{
			// Create new random elements
			if (_size == 0)
				elt = rand();
			else
				elt = static_cast<long>((double(rand())/RAND_MAX)*double(_size));

#ifdef TRACE
			double temp = elt;
			cout << "random double = " << temp 
			     << "    random NTL::RR = " << elt << endl;
#endif // TRACE

			return elt;
    
		} // element& operator() (void)



	//@} 
} // namespace LinBox

#endif // __LINBOX_field_ntl_rr_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
