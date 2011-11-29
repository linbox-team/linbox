
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
#include "linbox/linbox-config.h"
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

  /** 
  \brief give LinBox fields an allure of Givaro Fields
  \ingroup field

   *  This class adds the necessary requirements allowing 
   *  the construction of an extension of a LinBox field.
   */ 
    template< class BaseField >
    struct GivaroField : public BaseField 
    {
        typedef typename BaseField::Element TT;
        typedef typename Signed_Trait<TT>::unsigned_type UTT;
        typedef TT Rep;
        typedef GivaroField<BaseField> Self_t;
        typedef Rep Element;
        typedef UTT Residu_t;

        Element zero, one;
        GivaroField(const BaseField& bf) : BaseField(bf) {
            this->init(zero,0UL);
            this->init(one, 1UL);
        }


            // -- amxy: r <- c - a * b mod p
        Rep& amxy (Rep& r, const Rep a, const Rep b, const Rep c) const {
            Rep tmp;
            this->mul(tmp, a, b);
            this->assign(r,c);
            return this->subin(r,tmp);
        }

        bool areNEqual ( const Rep a, const Rep b) const {
            return ! this->areEqual(a,b);
        }

            // Access to the modulus, characteristic, size, exponent
        UTT residu() const { integer c; BaseField::characteristic(c); return UTT(c); }
        UTT characteristic() const  { integer c; BaseField::characteristic(c); return UTT(c); }
        UTT cardinality() const  { integer c; BaseField::cardinality(c); return UTT(c); }
        UTT exponent() const { return 1; }
        UTT size() const  { integer c; BaseField::cardinality(c); return UTT(c); }


            // ----- random generators
        template<class RandIter> Rep& random(RandIter& g, Rep& r) const { return r = g() ; }
        template<class RandIter> Rep& random(RandIter& g, Rep& r, long s) const { return r = g() ; }
        template<class RandIter> Rep& random(RandIter& g, Rep& r, const Rep& b) const { return r = g() ; }
        template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r) const { return r = g() ; }
        template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r, long s) const { return r = g() ; }
       template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r, const Rep& b) const { return r = g() ; }

    };
        
        

  /** This template class is define just to be in phase with the LinBox
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
    template<>
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
#endif // __FIELD_GIVARO_Extension
