/* File: src/library/objects/field/param_modular.h
 * Author: William Turner for the LinBox group
 */

#ifndef _PARAM_MODULAR_
#define _PARAM_MODULAR_

#include <iostream>
#include "LinBox/integer.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 

  // Forward declarations
  class param_modular_randIter;

  /** Parameterized field modulo prime number.
   * The primality of the modulus will not be checked, so 
   * it is the programmer's responsibility to supply a prime modulus.
   * This class implements 
   * a field of unparamterized integers modulo a prime integer.
   * Field has (non-static) member to contain modulus of field.
   */
  class param_modular
  {
  public:

    /** Element type.
     * It must meet the common object interface of elements as given in the
     * the archetype Element_archetype.
     */
    typedef integer element;

    /** Random iterator generator type.
     * It must meet the common object interface of random element generators
     * as given in the the archetype RandIter_archetype.
     */
    typedef param_modular_randIter randIter;

    /** @name Object Management
     */
    //@{
 
    /** Default constructor.
     */
    param_modular(void) {}

    /** Constructor from an integer.
     * Sets the modulus of the field throug the static member of the 
     * element type.
     * @param value constant reference to integer prime modulus
     */
    param_modular(const integer& value) : _modulus(value) {}

    /** Copy constructor.
     * Constructs param_modular object by copying the field.
     * This is required to allow field objects to be passed by value
     * into functions.
     * @param  F param_modular object.
     */
    param_modular(const param_modular& F) 
      : _modulus(F._modulus) {}
 
    /** Assignment operator.
     * Required by abstract base class.
     * @return reference to param_modular object for self
     * @param F constant reference to param_modular object
     */
    param_modular& operator= (const param_modular& F)
    { return *this; }

    /** Initialization of field base element from an integer.
     * Behaves like C++ allocator construct.
     * This function assumes the output field base element x has already been
     * constructed, but that it is not already initialized.
     * This is not a specialization of the template function because
     * such a specialization is not allowed inside the class declaration.
     * @return reference to field base element.
     * @param x field base element to contain output (reference returned).
     * @param y integer.
     */
    element& init(element& x, const integer& y) const
    { 
      x = y % _modulus;
      if (x < 0) x += _modulus;
      return x;
    }
 
   /** Conversion of field base element to a template class T.
     * This function assumes the output field base element x has already been
     * constructed, but that it is not already initialized.
     * @return reference to template class T.
     * @param x template class T to contain output (reference returned).
     * @param y constant field base element.
     */
    integer& convert(integer& x, const element& y) const
    { return x = y; }
 
    /** Assignment of one field base element to another.
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    element& assign(element& x, const element& y) const { return x = y; }

    /** Cardinality.
     * Return integer representing cardinality of the domain.
     * Returns a non-negative integer for all domains with finite
     * cardinality, and returns -1 to signify a domain of infinite
     * cardinality.
     * @return integer representing cardinality of the domain
     */
    const integer& cardinality(void) const
    { return *(new integer(_modulus)); }
 
    /** Characteristic.
     * Return integer representing characteristic of the domain.
     * Returns a positive integer to all domains with finite characteristic,
     * and returns 0 to signify a domain of infinite characteristic.
     * @return integer representing characteristic of the domain.
     */
    const integer& characteristic(void) const
    { return *(new integer(_modulus)); }

    //@} Object Management

    /** @name Arithmetic Operations
     * x <- y op z; x <- op y
     * These operations require all elements, including x, to be initialized
     * before the operation is called.  Uninitialized field base elements will
     * give undefined results.
     */
    //@{

    /** Equality of two elements.
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return boolean true if equal, false if not.
     * @param  x field base element
     * @param  y field base element
     */
    bool areEqual(const element& x, const element& y) const
    { return x == y; }

    /** Addition.
     * x = y + z
     * This function assumes all the field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     * @param  z field base element.
     */
    element& add(element& x, const element& y, const element& z) const
    { return x = (y + z) % _modulus; }
 
    /** Subtraction.
     * x = y - z
     * This function assumes all the field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     * @param  z field base element.
     */
    element& sub(element& x, const element& y, const element& z) const
    { 
      x = (y - z) % _modulus;
      if (x < 0) x += _modulus;
      return x;
    }
 
    /** Multiplication.
     * x = y * z
     * This function assumes all the field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     * @param  z field base element.
     */
    element& mul(element& x, const element& y, const element& z) const
    { return x = (y * z) % _modulus; }
 
    /** Division.
     * x = y / z
     * This function assumes all the field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     * @param  z field base element.
     */
    element& div(element& x, const element& y, const element& z) const
    { 
      element temp;
      inv(temp, z);
      return mul(x, y, temp);
    }
 
    /** Additive Inverse (Negation).
     * x = - y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    element& neg(element& x, const element& y) const
    { return x = _modulus - y; }
 
    /** Multiplicative Inverse.
     * x = 1 / y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    element& inv(element& x, const element& y) const
    {
      // The extended Euclidean algoritm
      integer x_int, y_int, q, tx, ty, temp;
      x_int = _modulus; 
      y_int = y;
      tx = 0; 
      ty = 1;

      while(y_int != 0)
      {
	// always: gcd(modulus,residue) = gcd(x_int,y_int)
        //         sx*modulus + tx*residue = x_int
        //         sy*modulus + ty*residue = y_int
        q = x_int / y_int; // integer quotient
        temp = y_int;  y_int  = x_int  - q*y_int;  x_int  = temp;
        temp = ty; ty = tx - q*ty; tx = temp;
      }

      // now x_int = gcd(modulus,residue)
      x = tx;
      if (x < 0) x += _modulus;

      return x;
    } // element& inv(element&, const element&) const

    //@} Arithmetic Operations
 
    /** @name Inplace Arithmetic Operations
     * x <- x op y; x <- op x
     */
    //@{

    /** Zero equality.
     * Test if field base element is equal to zero.
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return boolean true if equals zero, false if not.
     * @param  x field base element.
     */
    bool isZero(const element& x) const
    { return x == 0; }
 
    /** One equality.
     * Test if field base element is equal to one.
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return boolean true if equals one, false if not.
     * @param  x field base element.
     */
    bool isOne(const element& x) const
    { return x == 1; }

    /** Inplace Addition.
     * x += y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    element& addin(element& x, const element& y) const
    { 
      x += y;
      x %= _modulus;
      return x;
    }
 
    /** Inplace Subtraction.
     * x -= y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    element& subin(element& x, 
			    const element& y) const
    {
      x -= y;
      x %= _modulus;
      if (x < 0) x += _modulus;
      return x;
    }
 
    /** Inplace Multiplication.
     * x *= y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    element& mulin(element& x, 
			    const element& y) const
    {
      x *= y;
      x %= _modulus;
      return x;
    }
 
    /** Inplace Division.
     * x /= y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    element& divin(element& x, 
			    const element& y) const
    {
      element temp;
      inv(temp, y);
      x *= temp;
      x %= _modulus;
      return x;
    }
 
    /** Inplace Additive Inverse (Inplace Negation).
     * x = - x
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     */
    element& negin(element& x) const
    {
      x = _modulus - x;
      return x;
    }
 
    /** Inplace Multiplicative Inverse.
     * x = 1 / x
     * This function assumes the field base elementhas already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     */
    element& invin(element& x) const
    { return inv(x, x); }

    //@} Inplace Arithmetic Operations

    /** @name Input/Output Operations */
    //@{

    /** Print field.
     * @return output stream to which field is written.
     * @param  os  output stream to which field is written.
     */
    ostream& write(ostream& os) const 
    { return os << " mod " << _modulus; }
 
    /** Read field.
     * @return input stream from which field is read.
     * @param  is  input stream from which field is read.
     */
    istream& read(istream& is) { return is >> _modulus; }

    /** Print field base element.
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return output stream to which field base element is written.
     * @param  os  output stream to which field base element is written.
     * @param  x   field base element.
     */
    ostream& write(ostream& os, const element& x) const
    { return os << x << " mod " << _modulus; }
 
    /** Read field base element.
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return input stream from which field base element is read.
     * @param  is  input stream from which field base element is read.
     * @param  x   field base element.
     */
    istream& read(istream& is, element& x) const
    { 
      is >> x;

      x %= _modulus;
      if (x < 0) x += _modulus;

      return is; 
    }

    //@}

  private:

    /// Private (non-static) integer for modulus
    integer _modulus;

  }; // class param_modular

} // namespace LinBox

#include "LinBox/param_modular_randiter.h"

#endif // _PARAM_MODULAR_
