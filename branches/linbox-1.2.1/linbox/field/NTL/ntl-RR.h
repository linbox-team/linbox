/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/field/ntl-RR.h
 * Copyright (C) 1999-2005 William J Turner,
 *               2001 Bradford Hovinen
 * Copyright (C) 2011 LinBox
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

/*! @file field/NTL/ntl-RR.h
 * @ingroup field
 * @ingroup NTL
 * @brief NO DOC
 */

#ifndef __LINBOX_field_ntl_rr_H
#define __LINBOX_field_ntl_rr_H

#include "linbox/linbox-config.h"

#ifndef __LINBOX_HAVE_NTL
#error "you need NTL here"
#endif


#include <NTL/tools.h>
#include <NTL/RR.h>

#include "linbox/util/debug.h"

#include "linbox/field/unparametric.h"
#include "linbox/randiter/unparametric.h"
#include "linbox/field/field-traits.h"


// Namespace in which all LinBox library code resides
namespace LinBox
{

	class NTL_RR_Initialiser {
	public :
		NTL_RR_Initialiser () { }
	};



	/** @name class RR.
	  \brief
	 * Rational number field.
	 * This field is provided as a convenience in a few places.
	 * Use with caution because expression swell.
	 *
	 * This specialization allows the \ref UnparametricField template class to be
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
	NTL::RR& Caster(NTL::RR& x, const integer& y)
	{
		return x = NTL::to_RR(static_cast<const long&>(y));
	}
	template <>
	NTL::RR& Caster(NTL::RR& x, const double& y)
	{
		return x = NTL::to_RR((long)(y));
	}
	template <>
	NTL::RR& Caster(NTL::RR& x, const int& y)
	{
		return x = NTL::to_RR((long)(y));
	}

	template <>
	NTL::RR& Caster(NTL::RR& x, const long int& y)
	{
		return x = NTL::to_RR((long)(y));
	}



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
	integer& Caster(integer& x, const NTL::RR& y)
	{
		return x = static_cast<integer>(to_long(y));
	}

	struct NTL_RR: public NTL_RR_Initialiser, public FFPACK::UnparametricOperations<NTL::RR> {
		typedef NTL::RR Element ;
		typedef FFPACK::UnparametricOperations<Element> Father_t ;
		typedef UnparametricRandIter<Element> RandIter;

		const Element zero,one,mone ;

		NTL_RR() :
			NTL_RR_Initialiser(),Father_t ()
			,zero( NTL::to_RR(0)),one( NTL::to_RR(1)),mone(-one)
		{
			// no default - allow initialization of ZZ_p directly by user.
		}


		/** Multiplicative Inverse.
		 * x = 1 / y
		 * This function assumes both field elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element& inv(Element& x, const Element& y) const
		{
			return x = NTL::inv(y);
		}

		/** Zero equality.
		 * Test if field element is equal to zero.
		 * This function assumes the field element has already been
		 * constructed and initialized.
		 * In this specialization, NTL's IsZero function is called.
		 * @return boolean true if equals zero, false if not.
		 * @param  x field element.
		 */
		bool isZero(const Element& x) const
		{
			return static_cast<bool>(IsZero(x));
		}

		/** One equality.
		 * Test if field element is equal to one.
		 * This function assumes the field element has already been
		 * constructed and initialized.
		 * In this specialization, NTL's IsOne function is called.
		 * @return boolean true if equals one, false if not.
		 * @param  x field element.
		 */
		bool isOne(const Element& x) const
		{
			return static_cast<bool>(IsOne(x));
		}

		/** Inplace Multiplicative Inverse.
		 * x = 1 / x
		 * This function assumes both field elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 */
		Element& invin(Element& x) const
		{
			return x = NTL::inv(x);
		}

		/** Print field.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 */
		std::ostream& write(std::ostream& os) const
		{
			return os << "unparameterized field Element";
		}

		integer & cardinality(integer &c) const
		{
			return c = -1L;
		}

		integer & characteristic(integer &c) const
		{
			return c = 0UL;
		}

		std::ostream &write (std::ostream &os, const Element &x) const {
			return FFPACK::UnparametricOperations<Element>::write(os,x);
		}

		template <typename Src>
		Element& init (Element& x, const Src& s) const
		{
			return Caster (x, s);
		}

		template <typename T>
		T& convert (T &x, const Element &y) const
		{
			return Caster (x,y);
		}


	};

	template <>
	class UnparametricRandIter<NTL::RR> {
		typedef NTL::RR Element ;
	protected:
		integer _size,_seed;
	public:

		UnparametricRandIter<NTL::RR> (const NTL_RR & F,
					       const integer& size = 0,
					       const integer& seed = 0) :
			_size(size), _seed(seed)
		{
			if (_seed == integer(0)) _seed = integer(time(NULL));

			// integer cardinality;
			// F.cardinality(cardinality);
			// if (_size > cardinality)
				// _size = 0;

#ifdef TRACE
			std::cout << "created random generator with size " << _size
			<< " and seed " << _seed << std::endl;
#endif // TRACE

			// Seed random number generator
			NTL::SetSeed(NTL::to_ZZ(static_cast<long>(_seed)));
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
		Element& random(Element &elt) const
		{
			// Create new random elements
			if (_size == 0)
				elt = rand();
			else
				elt = static_cast<double>((double(rand())/RAND_MAX)*double(_size));

#ifdef TRACE
			double temp = elt;
			std::cout << "random double = " << temp << "    random Element = " << elt << std::endl;
#endif // TRACE

			return elt;

		} // element& operator() (void)

	};

	template <class Ring>
	struct ClassifyRing;

	template<class Element>
	struct ClassifyRing<UnparametricField<Element> >;

	template<>
	struct ClassifyRing<NTL_RR >{
		typedef RingCategories::ModularTag categoryTag;
	};
	//@}
} // namespace LinBox

#endif // __LINBOX_field_ntl_rr_H

