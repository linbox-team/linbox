/* linbox/field/givaro-gfq.h
 * Copyright (C) 2004 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

/* WARNING this wrapper works only with an improved version of Givaro.
 * This version of givaro won't be available for public yet.
 * But it is available on my web page.
 * You can send me a mail to get it or for others details.
 */

#ifndef __LINBOX_field_givaro_montgomery_H
#define __LINBOX_field_givaro_montgomery_H


#include <linbox/integer.h>
#include <linbox/field/field-interface.h>
#include <linbox/util/debug.h>
#include "linbox/linbox-config.h"
#include <linbox/field/field-traits.h>

//------------------------------------
// Files of Givaro library


#include <givaro/givmontg32.h>
#include <givaro/giv_randiter.h>
//------------------------------------

// Namespace in which all LinBox code resides
namespace LinBox 
{ 


	template <class Ring>
	struct ClassifyRing;

	class GivaroMontg;

	template<>
	struct ClassifyRing<GivaroMontg> {
		typedef RingCategories::ModularTag categoryTag;
	};
	/** 
	 \brief wrapper of Givaro's Montgomery<Std32>.
	 \ingroup field

	 *  This class is a modular representation with a Montgomery reduction
	 */   
	class GivaroMontg : public Montgomery<Std32>, public FieldInterface
	{
 
	public:

		/** Element type.
		 *  This type is inherited from the Givaro class Montgomery<Std32>
		 */
		typedef  Montgomery<Std32>::Rep Element;
    
		/** RandIter type
		 *  This type is inherited from the Givaro class Montgomery<Std32>
		 */	
		typedef GIV_randIter< Montgomery<Std32>, LinBox::integer >  RandIter;

		/** Constructor from an integer
		 *  this constructor use the ZpzDom<TAG> constructor
		 */
		GivaroMontg(const integer& p) :
			Montgomery<Std32>(static_cast<uint32>(long(p))) { }

		/** Constructor from an integer (takes degree of extension as 2nd parameter, must be 1)
		 *  this constructor use the ZpzDom<TAG> constructor
		 */
	  	GivaroMontg(const integer& p, const integer& k) :
			Montgomery<Std32>(static_cast<uint32>(long(p))) { 
			
			if (k!=1)
				throw PreconditionFailed(__FUNCTION__,__LINE__,"exponent must be 1");
		}
    
		/** Characteristic.
		 * Return integer representing characteristic of the domain.
		 * Returns a positive integer to all domains with finite characteristic,
		 * and returns 0 to signify a domain of infinite characteristic.
		 * @return integer representing characteristic of the domain.
		 */
		integer& characteristic(integer& c) const
		{return c=integer(static_cast<long>(Montgomery<Std32>::characteristic()));}

		long characteristic() const
		{return static_cast<long>(Montgomery<Std32>::characteristic());}
    
      
		/** Cardinality. 
		 * Return integer representing cardinality of the domain.
		 * Returns a non-negative integer for all domains with finite
		 * cardinality, and returns -1 to signify a domain of infinite
		 * cardinality.
		 * @return integer representing cardinality of the domain
		 */
		integer& cardinality(integer& c) const
		{ return c=integer(static_cast<long>(Montgomery<Std32>::size()));}
 

		/** Initialization of field base Element from an integer.
		 * Behaves like C++ allocator construct.
		 * This function assumes the output field base Element x has already been
		 * constructed, but that it is not already initialized.
		 * We assume that the type of Element is short int.
		 * this methos is just a simple cast.
		 * @return reference to field base Element.
		 * @param x field base Element to contain output (reference returned).
		 * @param y integer.
		 */  
		Element& init(Element& x , const integer& y=0) const
		{ 			
			return Montgomery<Std32>::init( x,long(y % (integer)_p));
		}
      
		Element& init(Element& x , const double y) const
		{ return Montgomery<Std32>::init( x, y);}

		/** Conversion of field base element to an integer.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.
		 * @return reference to an integer.
		 * @param x integer to contain output (reference returned).
		 * @param y constant field base element.
		 */
		integer& convert(integer& x, const Element& y) const
		{
			long tmp;	
			//	return x = *(new integer(Montgomery<Std32>::convert(tmp,y)));
			return x = integer(Montgomery<Std32>::convert(tmp,y));
		}
		double& convert(double& x, const Element& y) const
		{
			return Montgomery<Std32>::convert( x, y);
		}

		//bool isZero(const Element& x) const { return Montgomery<Std32>::isZero(x); }


		static inline int getMaxModulus() { return 40504; }

	}; // class GivaroMontg
 


} // namespace LinBox

#endif // __LINBOX_field_givaro_montgomery_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
