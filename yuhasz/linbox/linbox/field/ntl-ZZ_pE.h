#include <linbox/field/unparametric.h>
#include <linbox/randiter/unparametric.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ.h>
#include <time.h>
#include "linbox-config.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

using LinBox::Reader;
using LinBox::Writer;

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using std::istream;
using std::ostream;
using std::stringstream;
using std::ostringstream;
using std::string;
using std::vector;

#endif

namespace LinBox
{
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
	UnparametricRandIter(Reader &R) {
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
      Element& random (Element& x)
	{
	   NTL::random(x);
	   return x;
	}

#ifdef __LINBOX_XMLENABLED
      ostream &write(ostream &os) const
      {
	      Writer W;
	      if( toTag(W))
		      W.write(os);

	      return os;
      }

      bool toTag(Writer &W) const
      {
	      string s;
	      W.setTagName("randiter");
	      W.setAttribute("seed", Writer::numToString(s, _seed));
	      W.setAttribute("size", Writer::numToString(s, _size));

	      return true;
      }
#endif


    protected:
      size_t _size;
      size_t _seed;
    };
}

namespace LinBox
{


#ifdef __LINBOX_XMLENABLED
  template<>
    UnparametricField<NTL::ZZ_pE>::UnparametricField(Reader &R) 
    {
	    ostringstream oss;
	    string s;
	    size_t i;
	    long e;
	    NTL::ZZ m;
	    NTL::ZZ_pX poly;
	    vector<NTL::ZZ_p> v;

	    if(!R.expectTagName("field")) return;
	    if(!R.expectAttributeNum("cardinality", _card)) return;
	    if(!R.expectChildTag()) return;
	    R.traverseChild();

	    if(!R.expectTagName("finite") || !R.expectChildTag()) return;

	    R.traverseChild();
	    if(!R.expectTagName("characteristic") || !R.expectChildTag()) return;
	    R.traverseChild();
	    if(!R.expectTagNum(m)) return;
	    oss << m;
	    _p = Integer(oss.str().c_str());

	    NTL::ZZ_p::init(m);

	    R.upToParent();
	    R.upToParent();

	    if(!R.getNextChild()) {
		    R.setErrorString("finite field did not have extension or polynomial modulus.");
		    R.setErrorCode(Reader::OTHER);
		    return;
	    }
	    R.traverseChild();
	    if(!R.expectTagName("extension") || !R.expectChildTag()) return;
	    R.traverseChild();
	    if(!R.expectTagNum(e)) return;
	    R.upToParent();
	    R.upToParent();

	    if(!R.getNextChild()) {
		    // usually we'd insert code here for generating a default
		    // polynomial to act as a modulus, however we have none here,
		    // so instead we simply return an error
		    R.setErrorString("Got finite field with characteristic and extension degree, but no polynomial modulus.");
		    R.setErrorCode(Reader::OTHER);
		    return;
	    }
	    R.traverseChild();
	    if(!R.expectTagName("polynomial") || !R.expectTagNumVector(v)) return;
	    vector<NTL::ZZ_p>::const_iterator iter = v.begin();
	    i = 0;
	    while(iter != v.end()) {
		    
		    SetCoeff(poly, i, *iter);
		    ++i; ++iter;
	    }

	    // finally, initalize ZZ_pE
	    NTL::ZZ_pE::init(poly);

	    return;
    }

#endif




  template<>
    NTL::ZZ_pE& UnparametricField<NTL::ZZ_pE>::init (NTL::ZZ_pE &x, const integer &y) const
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
  NTL::ZZ_pE& UnparametricField<NTL::ZZ_pE>::inv(NTL::ZZ_pE& x, const NTL::ZZ_pE& y) const
    {
      x=NTL::to_ZZ_pE(1)/y;
      return x;
    }
   NTL::ZZ_pE& UnparametricField<NTL::ZZ_pE>::invin(NTL::ZZ_pE& x) const
     {
       x=NTL::to_ZZ_pE(1)/x;
       return x;
     }

#ifdef __LINBOX_XMLENABLED

   template <>
   bool UnparametricField<NTL::ZZ_pE>::toTag(Writer &W) const
   {
	   string s;
	   NTL::ZZ_pX poly = NTL::ZZ_pE::modulus();
	   long i;

	   W.setTagName("field");
	   W.setAttribute("implDetail", "ntl-ZZpE");
	   W.setAttribute("cardinality", Writer::numToString(s, _card));

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
   ostream &UnparametricField<NTL::ZZ_pE>::write(ostream &os) const
   {
	   Writer W;
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
   bool UnparametricField<NTL::ZZ_pE>::toTag(Writer &W, const Element &e) const
   {
	   NTL::ZZ_pX poly = rep(e);
	   NTL::ZZ accum, base = NTL::ZZ_p::modulus();
	   long i;
	   string s;

	   accum = 0;
	   for(i = deg(poly); i >= 0; --i) {
		   accum *= base;
		   accum += rep(coeff(poly, i));
	   }


	   W.setTagName("cn");
	   W.addDataChild(Writer::numToString(s, accum));

	   return true;
   }

   template <>
   ostream &UnparametricField<NTL::ZZ_pE>::write(ostream &os, const Element &e) const
   {

	   Writer W;
	   if( toTag(W, e))
		   W.write(os);

	   return os;
   }



   template <>
   bool UnparametricField<NTL::ZZ_pE>::fromTag(Reader &R, Element &e) const
   {
	   NTL::ZZ total, base = NTL::ZZ_p::modulus(), rem;
	   stringstream ss;

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
   istream &UnparametricField<NTL::ZZ_pE>::read(istream &is, Element &e) const
   {
	   Reader R(is);
	   if( !fromTag(R, e)) {
		   is.setstate(istream::failbit);
		   if(!R.initalized()) {
			   is.setstate(istream::badbit);
		   }
	   }

	   return is;
   }


#endif
	   

		  
	   


}
