
/* -*- mode: c; style: linux -*- */

/* linbox/field/givaro-gfq.h
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 * JGD 12.06.2002 : -- I don't see the need of *(new in convert
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


#include "linbox/integer.h"

//------------------------------------
// Files of Givaro library
#ifdef __TIMER_H
#define _TIMER_H_
#else
#define __TIMER_H
#endif

#include <givgfq.h>
#include <giv_randiter.h>
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
 class GivaroGfq : public GFqDom<long>
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
    GivaroGfq(const integer& p, const integer& k) :
      GFqDom<long>(static_cast<UTT>(long(p)), static_cast<UTT>(long(k))) {}
    

    /** Characteristic.
     * Return integer representing characteristic of the domain.
     * Returns a positive integer to all domains with finite characteristic,
     * and returns 0 to signify a domain of infinite characteristic.
     * @return integer representing characteristic of the domain.
     */
    integer& characteristic(integer& c) const
      {return c=integer(static_cast<long>(GFqDom<long>::characteristic()));}
      
      
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

  }; // class GivaroGfq
 

} // namespace LinBox

#endif // __FIELD_GIVARO_GFQ
