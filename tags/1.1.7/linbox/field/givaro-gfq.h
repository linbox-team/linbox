/* linbox/field/givaro-gfq.h
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 * JGD 12.06.2002 : -- I don't see the need of *(new in convert
 * JGD 19.09.2003 : added isZero
 * WJT 24.06.2005 : Removed using declarations
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

#ifndef __LINBOX_field_givaro_gfq_H
#define __LINBOX_field_givaro_gfq_H


#include <linbox/integer.h>
#include <linbox/field/field-traits.h>
#include <linbox/field/field-interface.h>
#include <linbox/util/debug.h>
#include <linbox/field/hom.h>
#include "linbox/linbox-config.h"


//------------------------------------
// Files of Givaro library


#include <givaro/givgfq.h>
#include <givaro/giv_randiter.h>
#include <givaro/givpoly1factor.h>
//------------------------------------

// Namespace in which all LinBox code resides
namespace LinBox 
{ 

	template <class Ring>
	struct ClassifyRing;

	class GivaroGfq;

	template<>
	struct ClassifyRing<GivaroGfq> {
		typedef RingCategories::ModularTag categoryTag;
        };


	class GivaroGfq;
	
	template<>
	inline integer& FieldTraits<GivaroGfq>::maxModulus( integer& i )
		{ return i = integer( 32749 ); } // prevprime( 2^15 )

	template<>
	inline bool FieldTraits<GivaroGfq>::goodModulus( const integer& i ) {
		integer max;
		if( i < 2 || i > FieldTraits<GivaroGfq>::maxModulus(max) )
			return false;
		return probab_prime( i, 10 );
	}

	template<>
	inline integer& FieldTraits<GivaroGfq>::maxExponent( integer& i )
		{ return i = 20; } // Cardinality must be <= 2^20


  /** wrapper of Givaro's GFqDom<int32>  class 
  \ingroup field

   *  This class allows to construct only extension fields with a prime characteristic.
   */   
 class GivaroGfq : public GFqDom<int32>, public FieldInterface
  {
 
  public:

    /** Element type.
     *  This type is inherited from the Givaro class GFqDom<int32>
     */
    typedef  GFqDom<int32>::Rep Element;
    
    /** RandIter type
     *  This type is inherited from the Givaro class GFqDom<TAG>
     */	
    typedef GIV_randIter< GFqDom<int32>, LinBox::integer >  RandIter;

    /** Empty Constructor 
     */
    GivaroGfq() : GFqDom<int32>() { }

    /** Constructor from an integer
     *  this constructor use the ZpzDom<TAG> constructor
     */
    GivaroGfq(const integer& p, const integer& k=1) :
      GFqDom<int32>(static_cast<UTT>(int32(p)), static_cast<UTT>(int32(k))) {
	//enforce that the cardinality must be <2^16, for givaro-gfq
	int32 pl=p;
	for(int32 i=1;i<k;++i) pl*=(int32)p;
	if(!FieldTraits<GivaroGfq>::goodModulus(p)) 
		throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus be between 2 and 2^15 and prime");
	else if(pl>(1<<20)) throw PreconditionFailed(__FUNCTION__,__LINE__,"cardinality must be < 2^20");

	}

    // Dan Roche 6-15-04
    // This constructor takes a vector of ints that represent the polynomial
    // to use (for modular arithmetic on the extension field).
    // Mostly copied from givaro/givgfq.inl
    GivaroGfq(const integer& p, const integer& k, const std::vector<integer>& modPoly)
      : GFqDom<int32>(static_cast<UTT>(int32(p)), static_cast<UTT>(int32(k))) {

        //enforce that the cardinality must be <2^16, for givaro-gfq
        int32 pl=p;
        for(int32 i=1;i<k;++i) pl*=(int32)p;
        if(!FieldTraits<GivaroGfq>::goodModulus(p)) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus be between 2 and 2^15 and prime");
        else if(pl>=(1<<16)) throw PreconditionFailed(__FUNCTION__,__LINE__,"cardinality must be < 2^16");
	
	if( k < 2 ) throw PreconditionFailed(__FUNCTION__,__LINE__,"exponent must be >1 if polynomial is specified");

	if(modPoly.size() != (size_t)(k+1)) throw PreconditionFailed(__FUNCTION__,__LINE__,"Polynomial must be of order k+1");

	GFqDom<int32> Zp(p,1);
	typedef Poly1FactorDom< GFqDom<int32>, Dense > PolDom;
	PolDom Pdom( Zp );
	PolDom::Element Ft, F, G, H;

	Poly1Dom< GFqDom<int32>, Dense >::Rep tempVector(k+1);
	for( int i = 0; i < k+1; i++ )
		tempVector[i] = modPoly[i] % p;
	Pdom.assign( F, tempVector );

	Pdom.give_prim_root(G,F);
	Pdom.assign(H,G);

	typedef Poly1PadicDom< GFqDom<int32>, Dense > PadicDom;
	PadicDom PAD(Pdom);

	PAD.eval(_log2pol[1],H);
	for (UTT i = 2; i < _qm1; ++i) {
		Pdom.mulin(H, G);
		Pdom.modin(H, F);
		PAD.eval(_log2pol[i], H);
	}

	for (UTT i = 0; i < _q; ++i)
		_pol2log[ _log2pol[i] ] = 1;
	
	UTT a,b,r,P=p;
	for (UTT i = 1; i < _q; ++i) {
		a = _log2pol[i];
		r = a & P;
		if (r == (P - 1))
			b = a - r;
		else
			b = a + 1;
		_plus1[i] = _pol2log[b] - _qm1;
	}

	_plus1[_qm1o2] = 0;

    }

    /** Characteristic.
     * Return integer representing characteristic of the domain.
     * Returns a positive integer to all domains with finite characteristic,
     * and returns 0 to signify a domain of infinite characteristic.
     * @return integer representing characteristic of the domain.
     */
    integer& characteristic(integer& c) const
      {return c=integer(static_cast<int32>(GFqDom<int32>::characteristic()));}

    int32 characteristic() const
      {return static_cast<int32>(GFqDom<int32>::characteristic());}
    
      
    /** Cardinality. 
     * Return integer representing cardinality of the domain.
     * Returns a non-negative integer for all domains with finite
     * cardinality, and returns -1 to signify a domain of infinite
     * cardinality.
     * @return integer representing cardinality of the domain
     */
    integer& cardinality(integer& c) const
      { return c=integer(static_cast<int32>(GFqDom<int32>::size()));}
 

    integer cardinality() const
      { return integer(static_cast<int32>(GFqDom<int32>::cardinality()));}
 

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
    Element& init(Element& x , const integer& y = 0) const
	  { return GFqDom<int32>::init( x, int32(y % (integer) _q));}
      
          // TO BE OPTIMIZED
    Element& init(Element& x , const float y) const
      { return GFqDom<int32>::init( x, (double)y);}

      template<class YYY>
    Element& init(Element& x , const YYY& y) const
      { return GFqDom<int32>::init( x, y);}

     /** Conversion of field base Element to an integer.
     * This function assumes the output field base Element x has already been
     * constructed, but that it is not already initialized.
     * @return reference to an integer.
     * @param x integer to contain output (reference returned).
     * @param y constant field base Element.
     */
    integer& convert(integer& x, const Element& y) const
      {
	int32 tmp;	
	return x = integer(GFqDom<int32>::convert(tmp,y));
      }
          // TO BE OPTIMIZED
    float& convert(float& x, const Element& y) const
      {
          double tmp;
          GFqDom<int32>::convert( tmp, y);
          return x = (float)tmp;
      }

      template<class XXX>
    XXX& convert(XXX& x, const Element& y) const
      {
	return GFqDom<int32>::convert( x, y);
      }

    //bool isZero(const Element& x) const { return GFqDom<int32>::isZero(x); }


  }; // class GivaroGfq
 




    template<> 
    class Hom <GivaroGfq,GivaroGfq>
    {   public:
        typedef GivaroGfq Source;
        typedef GivaroGfq Target;
      
        typedef Source::Element SrcElt;
        typedef Target::Element Elt;
        
            //Hom(){}
            /**
             * Construct a homomorphism from a specific source ring S and target 
             * field T with Hom(S, T).  
             * Specializations define all actual homomorphisms.
             */
        Hom(const Source& S, const Target& T) : _source(S), _target(T){ }
        
            /** 
             * image(t, s) implements the homomorphism, assigning the 
             * t the value of the image of s under the mapping.
             *
             * The default behaviour goes through integers.
             */
        Elt& image(Elt& t, const SrcElt& s) {
            return _target.init(t, _source.convert(tmp,s));
        }
        
            /** If possible, preimage(s,t) assigns a value to s such that 
             * the image of s is t.  Otherwise behaviour is unspecified.
             * An error may be thrown, a conventional value may be set, or
             * an arb value set.
             *
             * The default behaviour goes through integers.
             */
        SrcElt& preimage(SrcElt& s, const Elt& t) {
            return _source.init(s, _target.convert(tmp,t));
        }
        const Source& source() { return _source;}
        const Target& target() { return _target;}
        
    private:
        integer tmp;
        Source _source;
        Target _target;
    }; // end Hom 



} // namespace LinBox

#endif // __LINBOX_field_givaro_gfq_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
