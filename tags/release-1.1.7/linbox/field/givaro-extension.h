/* linbox/field/givaro-gfq.h
 * Copyright (C) 2005 JGD
 *
 * Time-stamp: <22 Jun 10 10:02:34 Jean-Guillaume.Dumas@imag.fr>
 * ------------------------------------
 *
 * See COPYING for license information.
 */



#ifndef __LINBOX_field_givaro_extension_H
#define __LINBOX_field_givaro_extension_H


#include <linbox/integer.h>
#include <linbox/field/field-traits.h>
#include <linbox/field/field-interface.h>
#include <linbox/util/debug.h>
#include "linbox/linbox-config.h"
#include <linbox/field/field-traits.h>
#include <linbox/field/givaro-gfq.h>

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#include <iostream>
#include <string>
#include <vector>

#endif //__LINBOX_XMLENABLED

//---------------------------------------------
// Files of Givaro library
#include <givaro/givextension.h>
#include <givaro/giv_randiter.h>
//---------------------------------------------
// To convert linbox fields to Givaro interface
#include <linbox/field/givaro-field.h>

//---------------------------------------------
// Namespace in which all LinBox code resides
namespace LinBox 
{ 

	template <class Ring>
	struct ClassifyRing;

        template< class BaseField>
	class GivaroExtension;

#if !defined(__INTEL_COMPILER) && !defined(__CUDACC__)
	template<>
#endif
	template< class BaseField>
	struct ClassifyRing<GivaroExtension<BaseField> > {
		typedef RingCategories::ModularTag categoryTag;
        };

#if !defined(__INTEL_COMPILER) && !defined(__CUDACC__)
	template<>
#endif
    template< class BaseField>
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


  /** This template class is defined to be in phase with the LinBox
   *  archetype.
   *  Most of all methods are inherited from Extension  class
   *  of Givaro.
   *  these class allow to construct only extension field with a prime characteristic.
   */   
    template< class BaseField = GivaroGfq>
    class GivaroExtension : public Extension<GivaroField<BaseField> >, public FieldInterface
    {
 
        typedef GivaroExtension<GivaroField<BaseField> > Self_t;
  public:

    /** Element type.
     *  This type is inherited from the Givaro class Extension
     */
    typedef typename Extension<GivaroField<BaseField> >::Element Element;
    
    /** RandIter type
     *  This type is inherited from the Givaro class GFqDom<TAG>
     */	
    typedef GIV_ExtensionrandIter< Extension<GivaroField<BaseField> >, LinBox::integer >  RandIter;


    GivaroExtension() {}
        

    /** Constructor from an integer
     */
    GivaroExtension(const integer& p, const integer& k=1) :
      Extension<GivaroField<BaseField> >(static_cast<typename Extension<GivaroField<BaseField> >::Residu_t>(int32(p)), static_cast<typename Extension<GivaroField<BaseField> >::Residu_t>(int32(k))) {
    }

    /** Constructor extension of a base field 
     */
    GivaroExtension(const BaseField& bF, const integer& ext=1) :
      Extension<GivaroField<BaseField> >( GivaroField<BaseField>(bF), static_cast<typename Extension<GivaroField<BaseField> >::Residu_t>(int32(ext))) {
    }


    /** Copy Constructor 
     */
    GivaroExtension(const Self_t& F) :
      Extension<GivaroField<BaseField> >(F) {
    }

  }; // class GivaroExtension
 


  /** This template class is define just to be in phase with the LinBox
   *  archetype.
   *  Most of all methods are inherited from Extension  class
   *  of Givaro.
   *  these class allow to construct only extension field with a prime characteristic.
   */   
#ifndef __INTEL_COMPILER
    template<>
#endif
    class GivaroExtension<GivaroGfq> : public Extension<GFqDom<int32> >, public FieldInterface
    {
 
        typedef GivaroExtension<GivaroGfq> Self_t;
  public:

    /** Element type.
     *  This type is inherited from the Givaro class Extension
     */
    typedef Extension<GFqDom<int32> >::Element Element;
    
    /** RandIter type
     *  This type is inherited from the Givaro class GFqDom<TAG>
     */	
    typedef GIV_ExtensionrandIter< Extension< GFqDom<int32> >, LinBox::integer >  RandIter;

    /** Constructor from an integer
     */
    GivaroExtension(const integer& p, const integer& k=1) :
      Extension<GFqDom<int32> >(static_cast< Extension<GFqDom<int32> >::Residu_t>(int32(p)), static_cast< Extension<GFqDom<int32> >::Residu_t>(int32(k))) {
    }

    /** Constructor extension of a base field 
     */
    GivaroExtension(const GivaroGfq& bF, const integer& ext=1) :
      Extension<GFqDom<int32> >( static_cast< const Extension<GFqDom<int32> >::BaseField_t &>(bF), static_cast< Extension<GFqDom<int32> >::Residu_t>(int32(ext))) {
    }


    /** Copy Constructor 
     */
    GivaroExtension(const Self_t& F) :
      Extension<GFqDom<int32> >(F) {
    }

  }; // class GivaroExtension
 


} // namespace LinBox



// Specialization of homomorphism for basefield
#include "linbox/field/hom.h"
namespace LinBox 
{
    template< class BaseField>
    class Hom < BaseField, GivaroExtension<BaseField> >
	{
            typedef BaseField Source;
            typedef GivaroExtension<BaseField> Target;
        public:
		typedef typename Source::Element SrcElt;
		typedef typename Target::Element Elt;

		//Hom(){}
		/**
		 * Construct a homomorphism from a specific source ring S and target 
		 * field T with Hom(S, T).  The default behaviour is error.  
		 * Specializations define all actual homomorphisms.
		 */
		Hom(const Source& S, const Target& T) : _source(S), _target(T){}

		/** 
		 * image(t, s) implements the homomorphism, assigning the 
		 * t the value of the image of s under the mapping.
		 *
		 * The default behaviour is a no-op.
		 */
		Elt& image(Elt& t, const SrcElt& s) const {return _target.assign(t, s);}

		/** If possible, preimage(s,t) assigns a value to s such that 
		 * the image of s is t.  Otherwise behaviour is unspecified.
		 * An error may be thrown, a conventional value may be set, or
		 * an arb value set.
		 *
		 * The default behaviour is a no-op.
		 */
		SrcElt& preimage(SrcElt& s, const Elt& t) const {
//                     return _target.getEntry(s, Degree(0), t);
                   return _target.convert(s, t);
                }

		const Source& source() const { return _source;}
		const Target& target() const { return _target;}

	private:
		Source _source;
		Target _target;
    }; // end Hom 
}
#endif // __LINBOX_field_givaro_extension_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
