/* File: src/wrapper/by_library/Givaro/Field/lin_zpz_giv.h
 * Author: Pascal Giorgi for the LinBox group
 */

#ifndef _LIN_ZPZ_GIV_
#define _LIN_ZPZ_GIV_

//-------------------------------------
// Files of LinBox library
#include "LinBox/integer.h"

//-------------------------------------
// Files of Givaro library
#include <givzpz16std.h>
#include <givzpz32std.h>
#include <givzpz16table1.h>
#include <giv_randiter.h>

//--------------------------------------


// Namespace in which all LinBox code resides
namespace LinBox 
{ 


  /** This template class is define just to be in phase with the LinBox
   *  archetype.
   *  Most of all methods are inherited from ZpzDom<Std16>, ZpzDom<Std32>
   *  and ZpzDom<log16> class of Givaro.
   *  these class allow to construct only finite field with a prime modulus.
   */   

template <class TAG> class givaro_zpz : public ZpzDom<TAG>
  {
  public:

    /** Element type.
     *  This type is inherited from the Givaro class ZpzDom<TAG>
     */
    typedef typename ZpzDom<TAG>::Rep element;

    /** Constructor from an integer
     *  this constructor use the ZpzDom<TAG> constructor
     */
    givaro_zpz(const integer& p) : ZpzDom<TAG>(static_cast<int>(p)) {}
    

    /** Characteristic.
     * Return integer representing characteristic of the domain.
     * Returns a positive integer to all domains with finite characteristic,
     * and returns 0 to signify a domain of infinite characteristic.
     * @return integer representing characteristic of the domain.
     */
    integer& characteristic(integer& c) const
      {return c=*(new integer(static_cast<int>(ZpzDom<TAG>::size())));}
      
      
    /** Cardinality. 
     * Return integer representing cardinality of the domain.
     * Returns a non-negative integer for all domains with finite
     * cardinality, and returns -1 to signify a domain of infinite
     * cardinality.
     * @return integer representing cardinality of the domain
     */
    integer& cardinality(integer& c) const
      { return c=*(new integer(static_cast<int>(ZpzDom<TAG>::size())));}
 

   /** Conversion of field base element to an integer.
     * This function assumes the output field base element x has already been
     * constructed, but that it is not already initialized.
     * @return reference to an integer.
     * @param x integer to contain output (reference returned).
     * @param y constant field base element.
     */
    integer& convert(integer& x, const element& y) const
      {return x = integer(static_cast<int>(y));}
      

    /** Initialization of field base element from an integer.
     * Behaves like C++ allocator construct.
     * This function assumes the output field base element x has already been
     * constructed, but that it is not already initialized.
     * We assume that the type of element is short int.
     * this methos is just a simple cast.
     * @return reference to field base element.
     * @param x field base element to contain output (reference returned).
     * @param y integer.
     */  
    element& init(element& x , const integer& y=0) const
      { x= element(static_cast<int>(y)); return x;}
      
    
  }; // class givaro_zpz<TAG>
 

} // namespace LinBox

#endif // _LIN_ZPZ_GIV_
