
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/givaro-gfq.h
 * Copyright (C) 2005 JGD
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */



#ifndef __FIELD_GIVARO_EXTENSION
#define __FIELD_GIVARO_EXTENSION


#include <linbox/integer.h>
#include <linbox/field/field-traits.h>
#include <linbox/field/field-interface.h>
#include <linbox/util/debug.h>
#include "linbox-config.h"
#include <linbox/field/field-traits.h>
#include <linbox/field/givaro-gfq.h>

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#include <iostream>
#include <string>
#include <vector>

#endif

//------------------------------------
// Files of Givaro library


#include <givaro/givextension.h>
#include <givaro/giv_randiter.h>
//------------------------------------

// Namespace in which all LinBox code resides
namespace LinBox 
{ 

	template <class Ring>
	struct ClassifyRing;

        template< class BaseField>
	class GivaroExtension;

	template<> template< class BaseField>
	struct ClassifyRing<GivaroExtension<BaseField> > {
		typedef RingCategories::ModularTag categoryTag;
        };

    template<> template< class BaseField>
    struct FieldTraits< GivaroExtension<BaseField> >
    {
	typedef RingCategories::ModularTag categoryTag;
	
	static integer& maxModulus( integer& i )
            { return  FieldTraits<BaseField>::maxModulus(i); } 

	static bool goodModulus( const integer& i ) 
            { return  FieldTraits<BaseField>::goodModulus(i); } 
        
	// After that degree might not be correct ...
	static integer& maxExponent( integer& i ) { return i = 2147483648UL; }
	static bool goodExponent( const integer& i ) {
            integer max;
            return ( i >= 1 && i <= maxExponent( max ) );
	}
    };


  /** This template class is define just to be in phase with the LinBox
   *  archetype.
   *  Most of all methods are inherited from Extension  class
   *  of Givaro.
   *  these class allow to construct only extension field with a prime characteristic.
   */   
    template< class BaseField = GivaroGfq>
    class GivaroExtension : public Extension<BaseField>, public FieldInterface
    {
 
        typedef GivaroExtension<BaseField> Self_t;
  public:

    /** Element type.
     *  This type is inherited from the Givaro class Extension
     */
    typedef typename Extension<BaseField>::Element Element;
    
    /** RandIter type
     *  This type is inherited from the Givaro class GFqDom<TAG>
     */	
    typedef GIV_ExtensionrandIter< Extension<BaseField>, LinBox::integer >  RandIter;

    /** Constructor from an integer
     */
    GivaroExtension(const integer& p, const integer& k=1) :
      Extension<BaseField>(static_cast<typename Extension<BaseField>::Residu_t>(int32(p)), static_cast<typename Extension<BaseField>::Residu_t>(int32(k))) {
    }

    /** Constructor extension of a base field 
     */
    GivaroExtension(const BaseField& bF, const integer& ext=1) :
      Extension<BaseField>(bF, static_cast<typename Extension<BaseField>::Residu_t>(int32(ext))) {
    }


    /** Copy Constructor 
     */
    GivaroExtension(const Self_t& F) :
      Extension<BaseField>(F) {
    }

  }; // class GivaroExtension
 


} // namespace LinBox

#endif // __FIELD_GIVARO_Extension
