/* File: /src/wrappers/by_scope/field/ALP_param_field.h
 * Author: Pascal Giorgi for the Linbox group
 */

#ifndef  __ALP_FIELD
#define  __ALP_FIELD

//-----------------------------------
// Files of C/C++ library
#include "sstream.h"
#include <iostream.h>

//-----------------------------------
// Files of ALP library
#include <arithm/Zp.H>
#include <bignum.H>

//-----------------------------------
// Files of LinBox library
#include "LinBox/integer.h"
#include "unparam_field.h"

//-----------------------------------

// Namespace in which all LinBox library code resides
namespace LinBox
{


  /** @name ALP Fields. We assume that the ALP field is a prime field.
   * instanciations for wrapping object from 
   * \URL[Bernard Mourrain]{http://www-sop.inria.fr/saga/logiciels/ALP/}'s
   * objects to comply with the common object interface for 
   * \Ref{Linbox} fields.
   */


  /** The parameter (int p) defined in the templated class below is the modulus
   *  of the ALP_field
   *  Each methods are inherited from the unparam_field wrapper except
   *  default constructor, convert, cardinality, characteristic.
   */
   
  template <int p> class alp_zpz : public unparam_field<Z<p> > 
    {
    public:
      
       /** Element type.
       * It must meet the common object interface of elements as given in the 
       * archetype Element_archetype.
       */
      typedef Z<p> element;

      /** Random element generator.
       */
      typedef unparam_field<Z<p> >::randIter randIter;
      

      /** Default constructor
       */
      alp_zpz(void) { }



      /** Conversion of field base element to a template class T.
      * This function assumes the output field base element x has already been
      * constructed, but that it is not already initialized.
      * @return reference to integer.
      * @param x template class T to contain output (reference returned).
      * @param y constant field base element.
      */
      integer& convert(integer& x, const element& y) const
       {
	 
	 ostringstream o;
	 o<<y;
	 istringstream i(o.str());
	 char s[1024];
	 i>>s; 
	 return x = integer(s);
       }


      /** Cardinality.
       * Return integer representing cardinality of the domain.
       * Returns a non-negative integer for all domains with finite
       * cardinality, and returns -1 to signify a domain of infinite
       * cardinality.
       * @return integer representing cardinality of the domain
       */
      integer& cardinality(integer& c) const
	{ return c=*(new integer(p)); }
      
      
      /** Characteristic.
       * Return integer representing characteristic of the domain.
       * Returns a positive integer to all domains with finite characteristic,
       * and returns 0 to signify a domain of infinite characteristic.
       * @return integer representing characteristic of the domain.
       */
      integer& characteristic(integer& c) const
	{ return c=*(new integer(p)); }
      
    };// template<int p> class alp_zpz
}//  namespace Linbox
      
   
#endif // ALP_FIELD
