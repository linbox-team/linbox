/* Copyright (C) LinBox
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_field_ntl_zz_pe_H
#define __LINBOX_field_ntl_zz_pe_H

#include <linbox/field/unparametric.h>
#include <linbox/randiter/unparametric.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ.h>
#include <time.h>
#include "linbox/linbox-config.h"
#include <linbox/field/field-traits.h>

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#endif //__LINBOX_XMLENABLED

namespace LinBox
{
	
	template <class Ring>
	struct ClassifyRing;

	template<>
	struct ClassifyRing<UnparametricRandIter<NTL::ZZ_pE> > {
		typedef RingCategories::ModularTag categoryTag;
	};

  template<>
  class UnparametricRandIter<NTL::ZZ_pE>
    {
    public:
      typedef NTL::ZZ_pE Element;
      UnparametricRandIter<NTL::ZZ_pE>(const UnparametricField<NTL::ZZ_pE>& F =UnparametricField<NTL::ZZ_pE>(), 
				       const size_t& size = 0,
				       const size_t& seed = 0
				       )
	: _size(size), _seed(seed)
	{
	  if(_seed == 0)
	    NTL::SetSeed(NTL::to_ZZ(time(0)));
	  else
	    NTL::SetSeed(NTL::to_ZZ(_seed));
	}

#ifdef __LINBOX_XMLENABLED
	UnparametricRandIter(LinBox::Reader &R) {
		if(!R.expectTagName("randiter")) return;
		if(!R.expectAttributeNum("seed", _seed) || !R.expectAttributeNum("size", _size)) return;

		if(_seed == 0) _seed = time(NULL);

		NTL::SetSeed(NTL::to_ZZ(_seed));
	}
#endif

      
      UnparametricRandIter<NTL::ZZ_pE>(const UnparametricRandIter<NTL::ZZ_pE>& R)
	: _size(R._size), _seed(R._seed) 
	
	{
	  if(_seed == 0)
	    NTL::SetSeed(NTL::to_ZZ(time(0)));
	  else
	    NTL::SetSeed(NTL::to_ZZ(_seed));
	}
      Element& random (Element& x) const 
	{
	   NTL::random(x);
	   return x;
	}

    protected:
      size_t _size;
      size_t _seed;
    };
}

namespace LinBox
{

  /*
   * Define a parameterized class to handle easily UnparametricField<NTL::ZZ_pE> field
   */
  class NTL_ZZ_pE : public UnparametricField<NTL::ZZ_pE>
    {
    public:
      NTL_ZZ_pE (const integer &p, const integer &k) {
	
	NTL::ZZ_p::init(NTL::to_ZZ(std::string(p).data()));
	NTL::ZZ_pX irredPoly = NTL::BuildIrred_ZZ_pX ((long) k);
	NTL::ZZ_pE::init(irredPoly);
      }
      
    }; // end o class NTL_ZZ_pE
  




  template<>
    NTL::ZZ_pE& UnparametricField<NTL::ZZ_pE>::init (NTL::ZZ_pE &x, const integer &y) const
    {
      x=NTL::to_ZZ_pE(static_cast<long>(y));
      return x;
    }
   template<>
    NTL::ZZ_pE& UnparametricField<NTL::ZZ_pE>::init (NTL::ZZ_pE &x, const double &y) const
    {
      x=NTL::to_ZZ_pE(static_cast<long>(y));
      return x;
    }
 
  template<>
    bool UnparametricField<NTL::ZZ_pE>::isZero (const NTL::ZZ_pE& a) const
    {
      return NTL::IsZero(a);
    }
  
  template<>
    bool UnparametricField<NTL::ZZ_pE>::isOne (const NTL::ZZ_pE& a) const
    {
      return NTL::IsOne(a);
    }

  // Rich Seagraves, 7-15-03
  // On the orders of Dr Saunders, I'm re-writing init & convert so that
  // they convert a ZZpE into a padic number, ie a0 + a1x + a2x^2 +... ->
  // a0 + a1*p + a2*p^2 + ...
  //
  template<>
    integer& UnparametricField<NTL::ZZ_pE>::convert(integer& c, const NTL::ZZ_pE& e) const
    {
	    NTL::ZZ_pX poly = rep(e);
	    Integer base = _p;
	    long i;

	    c = 0;
	    for(i = deg(poly); i >= 0; --i) {
		    c *= base;
		    c +=  NTL::to_long(rep(coeff(poly, i)));
	    }

	    return c;
    }
  
  template<>
    integer& UnparametricField<NTL::ZZ_pE>::characteristic (integer &c) const
    {
      return c=static_cast<integer>(to_long(NTL::ZZ_p::modulus()));
      //NTL::ZZ_p::modulus();
    }
  
  template<>
    integer& UnparametricField<NTL::ZZ_pE>::cardinality(integer& c) const
    {
      c=static_cast<integer>(to_long(NTL::ZZ_p::modulus()));
      c=pow(c,NTL::ZZ_pE::degree());
      return c;
    }
  template<>
    NTL::ZZ_pE& UnparametricField<NTL::ZZ_pE>::inv(NTL::ZZ_pE& x, const NTL::ZZ_pE& y) const
    {
      x=NTL::to_ZZ_pE(1)/y;
      return x;
    }
  template<>
    NTL::ZZ_pE& UnparametricField<NTL::ZZ_pE>::invin(NTL::ZZ_pE& x) const
     {
       x=NTL::to_ZZ_pE(1)/x;
       return x;
     }

   template<>
     std::istream& UnparametricField<NTL::ZZ_pE>::read(std::istream& is, NTL::ZZ_pE& x) const
     {
       long tmp;
       is>>tmp;
       x=NTL::to_ZZ_pE(tmp);
       return is;
     }		  
   
   

#ifdef __LINBOX_XMLENABLED

   template <>
   bool UnparametricField<NTL::ZZ_pE>::toTag(LinBox::Writer &W) const
   {
	   std::string s;
	   NTL::ZZ_pX poly = NTL::ZZ_pE::modulus();
	   long i;

	   W.setTagName("field");
	   W.setAttribute("implDetail", "ntl-ZZpE");
	   W.setAttribute("cardinality", LinBox::Writer::numToString(s, _card));

	   W.addTagChild();
	   W.setTagName("finite");

	   W.addTagChild();
	   W.setTagName("characteristic");
	   W.addNum(_p);
	   W.upToParent();

	   W.addTagChild();
	   W.setTagName("extension");
	   W.addNum(deg(poly) + 1);
	   W.upToParent();

	   W.addTagChild();
	   W.setTagName("polynomial");
	   
	   for(i = 0; i <= deg(poly); ++i) {
		   W.addNum(coeff(poly, i));
	   }
	   W.upToParent();
	   W.upToParent();
	   W.upToParent();

	   return true;
   }

   template <> 
   std::ostream &UnparametricField<NTL::ZZ_pE>::write(std::ostream &os) const
   {
	   LinBox::Writer W;
	   if( toTag(W) )
		   W.write(os);

	   return os;
   }


   // Elemnt Reading & writing functions
   // BIG NOTE:  It was decided that for extension fields, the elements
   // would be represented using a single number that has the following 
   // property:  for an element e in ZZp[x], with e = a0 + a1x + a2x^2 + ...,
   // represent e as "<cn>n</cn>" where n = a0 + a1 * p + a2 * p^2 + ...
   //

   template <>
   bool UnparametricField<NTL::ZZ_pE>::toTag(LinBox::Writer &W, const Element &e) const
   {
	   NTL::ZZ_pX poly = rep(e);
	   NTL::ZZ accum, base = NTL::ZZ_p::modulus();
	   long i;
	   std::string s;

	   accum = 0;
	   for(i = deg(poly); i >= 0; --i) {
		   accum *= base;
		   accum += rep(coeff(poly, i));
	   }


	   W.setTagName("cn");
	   W.addDataChild(LinBox::Writer::numToString(s, accum));

	   return true;
   }

   template <>
   std::ostream &UnparametricField<NTL::ZZ_pE>::write(std::ostream &os, const Element &e) const
   {

	   LinBox::Writer W;
	   if( toTag(W, e))
		   W.write(os);

	   return os;
   }



   template <>
   bool UnparametricField<NTL::ZZ_pE>::fromTag(LinBox::Reader &R, Element &e) const
   {
	   NTL::ZZ total, base = NTL::ZZ_p::modulus(), rem;
	   std::stringstream ss;

	   if(!R.expectTagName("cn") || !R.expectChildTextNum(total))
		   return false;

	   ss << "[";
	   while(total > 0) {
		   rem = total % base;
		   total /= base;
		   ss << rem;
		   if(total > 0) ss << " ";
	   }

	   ss << "]";
	   
	   ss >> e; // use the extraction stream operator

	   return true;
   }

   template <>
   std::istream &UnparametricField<NTL::ZZ_pE>::read(std::istream &is, Element &e) const
   {
	   LinBox::Reader R(is);
	   if( !fromTag(R, e)) {
		   is.setstate(std::istream::failbit);
		   if(!R.initalized()) {
			   is.setstate(std::istream::badbit);
		   }
	   }

	   return is;
   }


#endif
	   


   
   
}

#endif //__LINBOX_field_ntl_zz_pe_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
