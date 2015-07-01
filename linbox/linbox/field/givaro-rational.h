
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/givaro-rational.h
 * Copyright (C) 2004 Gilles Villard
 *
 * Created  Gilles Villard < Gilles.Villard@ens-lyon.fr> 
 * ven oct  8 13:30:05 CEST 2004
 */


#ifndef __GIVARO_RATIONAL_H
#define __GIVARO_RATIONAL_H


#include "linbox/integer.h"
#include "linbox/field/unparametric.h"
#include <linbox/field/field-traits.h>


//------------------------------------
// Files of Givaro library

#include "givaro/givrational.h"
//------------------------------------


namespace LinBox 
{ 

	template <class Ring>
    struct ClassifyRing;

	class GivaroRational; 
	
	template<> 
	struct ClassifyRing<GivaroRational> {
		typedef RingCategories::RationalTag categoryTag;
	}; 

 class GivaroRational : public UnparametricField<Rational>
  {
 
  public:

    /** Element type.
     *  
     */
    typedef  Rational Element;
    

    Element& init(Element& x , const integer& y = 0) const
	  { return x=Rational(y);}




  }; // class 
 


} // namespace LinBox

#endif // __GIVARO_RATIONAL_H
