
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/givaro-gfq.h
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 * JGD 12.06.2002 : -- I don't see the need of *(new in convert
 * JGD 19.09.2003 : added isZero
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

#ifndef __FIELD_GIVARO_GFQ
#define __FIELD_GIVARO_GFQ


#include <linbox/integer.h>
#include <linbox/field/field-interface.h>
#include <linbox/util/debug.h>
#include "linbox-config.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

using LinBox::Reader;
using LinBox::Writer;

#include <iostream>
#include <string>

using std::istream;
using std::ostream;
using std::string;

#endif

//------------------------------------
// Files of Givaro library


#include <givaro/givgfq.h>
#include <givaro/giv_randiter.h>
//------------------------------------

// Namespace in which all LinBox code resides
namespace LinBox 
{ 


  /** This template class is define just to be in phase with the LinBox
   *  archetype.
   *  Most of all methods are inherited from GFqDom<long>  class
   *  of Givaro.
   *  these class allow to construct only extension field with a prime characteristic.
   */   
 class GivaroGfq : public GFqDom<long>, public FieldInterface
  {
 
  public:

    /** Element type.
     *  This type is inherited from the Givaro class GFqDom<long>
     */
    typedef  GFqDom<long>::Rep Element;
    
    /** RandIter type
     *  This type is inherited from the Givaro class GFqDom<TAG>
     */	
    typedef GIV_randIter< GFqDom<long>, LinBox::integer >  RandIter;

    /** Constructor from an integer
     *  this constructor use the ZpzDom<TAG> constructor
     */
    GivaroGfq(const integer& p, const integer& k=1) :
      GFqDom<long>(static_cast<UTT>(long(p)), static_cast<UTT>(long(k))) {

	//enforce that the cardinality must be <2^16, for givaro-gfq
	long pl=p;
	// Rich Seagraves 7-16-03: Line removed to take care of compile warning
	//	long kl=k; 
	for(long i=1;i<k;++i) pl*=(long)p;
	if(p<=1) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus  must be >1");
	else if(pl>=(1<<16)) throw PreconditionFailed(__FUNCTION__,__LINE__,"cardinality must be < 2^16");

	}
    
#ifdef __LINBOX_XMLENABLED
    // XML Reader constructor
    GivaroGfq(Reader &R)
    {
	    integer p, n;

	    if(!R.expectTagName("field") || !R.expectChildTag()) return;
	    R.traverseChild();

	    if(!R.expectTagName("finite") || !R.expectChildTag()) return;
	    R.traverseChild();

	    if(!R.expectTagName("characteristic") || !R.expectChildTag()) return;
	    R.traverseChild();
	    if(!R.expectTagNum(p));
	    R.upToParent();

	    if(R.getNextChild()) {
		    if(!R.expectChildTag()) return;
		    R.traverseChild();

		    if(!R.expectTagName("extension") || !R.expectChildTag()) return;
		    R.traverseChild();
		    if(!R.expectTagNum(n)) return;
		    R.upToParent();

		    R.upToParent();
		    R.getPrevChild();
	    }
	    else {
		    n = Integer(1);
	    }
	    R.upToParent();
	    R.upToParent();

	    // now try building using the above constructor.  Note, NO
	    // ATTEMPT IS MADE TO CATCH THE ERROR THIS METHOD CAN THROW, 
	    // IT IS ALLOWED TO PASS THROUGH
	    //
	    GivaroGfq oth(p, n);
	    *this = oth;

	    return;
    }
#endif


    /** Characteristic.
     * Return integer representing characteristic of the domain.
     * Returns a positive integer to all domains with finite characteristic,
     * and returns 0 to signify a domain of infinite characteristic.
     * @return integer representing characteristic of the domain.
     */
    integer& characteristic(integer& c) const
      {return c=integer(static_cast<long>(GFqDom<long>::characteristic()));}

    long characteristic() const
      {return static_cast<long>(GFqDom<long>::characteristic());}
    
      
    /** Cardinality. 
     * Return integer representing cardinality of the domain.
     * Returns a non-negative integer for all domains with finite
     * cardinality, and returns -1 to signify a domain of infinite
     * cardinality.
     * @return integer representing cardinality of the domain
     */
    integer& cardinality(integer& c) const
      { return c=integer(static_cast<long>(GFqDom<long>::size()));}
 

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
      { return GFqDom<long>::init( x,long(y));}
      
    Element& init(Element& x , const double y=0.0) const
      { return GFqDom<long>::init( x, y);}

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
//	return x = *(new integer(GFqDom<long>::convert(tmp,y)));
	return x = integer(GFqDom<long>::convert(tmp,y));
      }
    double& convert(double& x, const Element& y) const
      {
	return GFqDom<long>::convert( x, y);
      }

    bool isZero(const Element& x) const { return GFqDom<long>::iszero(x); }


#ifdef __LINBOX_XMLENABLED

	  ostream &write(ostream &os) const
	  {
		  Writer W;
		  if( toTag(W) )
			  W.write(os);

		  return os;
	  }

	  bool toTag(Writer &W) const
	  {
		  string s;
		  long card = GFqDom<long>::size();
		  size_t i = 0;

		  W.setTagName("field");
		  W.setAttribute("implDetail", "givaro-gfq");
		  W.setAttribute("cardinality", Writer::numToString(s, card));

		  W.addTagChild();
		  W.setTagName("finite");

		  W.addTagChild();
		  W.setTagName("characteristic");
		  W.addNum(GFqDom<long>::characteristic());
		  W.upToParent();
		  W.addTagChild();
		  W.setTagName("extension");

		  while(card > 1) {
			  card /= GFqDom<long>::characteristic();
			  ++i;
		  }
		  W.addNum(i);
		  W.upToParent();

		  W.upToParent();

		  return true;
	  }


	  // Special Note:  In LinBox, all elements of a field will be written
	  // in the following manner:  for e in ZZp[x] with 
	  // e = a0 + a1x + a2x^2 + ..., e is represented as:
	  // "<cn>n</cn>" where n = a0 + a1 * p + a2 * p^2 + ...
	  // 
	  ostream &write(ostream &os, const Element &e) const
	  {
		  Writer W;
		  if( toTag(W, e))
			  W.write(os);

		  return os;
	  }

	  bool toTag(Writer &W, const Element &e) const
	  {
		  string s;
		  long rep = _log2pol[ (unsigned long) e];

		  W.setTagName("cn");
		  W.addDataChild(Writer::numToString(s, rep));
		  
		  return true;
	  }

	  istream &read(istream &is, Element &e) const
	  {
		  Reader R(is);
		  if( !fromTag(R, e)) {
			  is.setstate(istream::failbit);
			  if(!R.initalized())
				  is.setstate(istream::badbit);
		  }

		  return is;
	  }

	  bool fromTag(Reader &R, Element &e) const
	  {
		  unsigned long i;

		  if(!R.expectTagName("cn") || !R.expectChildTextNum(i))
			  return false;

		  e = _pol2log[i];
		  return true;
	  }
			  

#endif


  }; // class GivaroGfq
 


} // namespace LinBox

#endif // __FIELD_GIVARO_GFQ
