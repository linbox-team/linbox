/* File: src/library/objects/field/modular.h
 * Author: William J. Turner for the LinBox group
 */
 
#ifndef _MODULAR_FIELD_
#define _MODULAR_FIELD_

#include <iostream>
#include "LinBox/integer.h"
#include "LinBox/unparam_field.h"
#include "LinBox/unparam_randiter.h"
#include "LinBox/faxpy.h"

// Namespace in which all LinBox code resides
namespace LinBox
{

  /** Unparameterized field modulo prime number.
   * The primality of the modulus will not be checked, so 
   * it is the programmer's responsibility to supply a prime modulus.
   * Used with \Ref{unparam_field} to create implement a \Ref{LinBox}
   * unparamterized field modulo a prime number.
   * Field has static member to contain modulus of field.
   */
  class modular
  {
  public:

    /** @name Object Management */
    //@{

    /** Default constructor.
     * Sets modulus to zero.
     */
    modular() : residue(0) {}
 
    /** Constructor from integer.
     * Sets residue to representative of residue class of integer.
     * @param  start prime integer modulus
     */
    modular(long start)
    {
      residue = start % modulus;
      if (residue < 0)
  	residue = residue + modulus;
    }
 
    /** Copy constructor.
     * Copies residue.
     * @param  z modular object.
     */
    modular(const modular& z) { residue = z.residue; }

    /** Assignment Operator.
     * "Copies" into *this, which is an existing object.
     * If not explicitly defined, each member is assigned.
     * @return reference to self
     * @param  z modular object
     */
    modular& operator=(const modular& z)
    {
      if (this != &z ) residue = z.residue; // guard against self assignment
      return *this;
    }
 
    /** Retrieve residue class.
     * @return integer representing residue class
     */
    long get_residue() const { return residue; }

    /** Change residue class.
     * Changes residue class of modular object to that of integer.
     * @param  start integer to which to change the residue class
     */
    void put_residue(long start)
    {
      residue = start % modulus;
      if (residue < 0)
  	residue = residue + modulus;
    }
 
    /** Retreive prime modulus.
     * @return prime modulus
     */
    static long get_modulus() { return modulus; }
 
    /** Change prime modulus.
     * Changes modulus of field to prime integer.
     * Static modulus means only one field can be implemented at a time.
     * @param  start prime modulus
     */
    static void put_modulus(long start) { modulus = start; }

    //@} Object Management

    /** @name Arithmetic Operators */
    //@{

    /** Equality.
     * @return boolean true if elements are equal, false if not.
     * @param  z modular object.
     */
    bool operator==(const modular& z) const { return residue == z.residue; }

    /** Addition.
     * @return self + z.
     * @param  z modular object
     */
    modular& operator+(const modular& z) const
    {
      modular* x_ptr = new modular;
      x_ptr->residue = ( residue + z.residue ) % modulus;
      return *x_ptr;
    }

    /** Subtraction.
      * @return self - z.
      * @param  z modular object
      */
    modular& operator-(const modular& z) const
    {
      modular* x_ptr = new modular;
      x_ptr->residue = ( residue - z.residue ) % modulus;
 
      if (x_ptr->residue < 0)
  	x_ptr->residue = x_ptr->residue + modulus;

      return *x_ptr;
    }
 
    /** Multiplication.
      * @return self * z.
      * @param  z modular object
      */
    modular& operator*(const modular& z) const
    {
      modular* x_ptr = new modular;
      x_ptr->residue = ( residue * z.residue ) % modulus;
      return *x_ptr;
    }

    /** Division.
      * @return self / z.
      * @param  z modular object
      */
    modular& operator/(const modular& z) const { return (*this) * z.recip(); }

    /** Negation.
      * @return - self
      */
    modular& operator-(void) const
    {
      modular* x_ptr = new modular;
      x_ptr->residue = modulus - residue;
      return *x_ptr;
    }

    //@} Arithmetic Operators

    /** @name Inplace Arithmetic Operators */
    //@{

    /** Inplace Addition.
      * self = self + z.
      * @return reference to self.
      * @param  z modular object
      */
    modular& operator+=(const modular& z) { return (*this) = (*this) + z; }

    /** Inplace Subtraction.
      * self = self - z.
      * @return reference to self.
      * @param  z modular object
      */
    modular& operator-=(const modular& z) { return (*this) = (*this) - z; }

    /** Inplace Multiplication.
      * self = self * z.
      * @return reference to self.
      * @param  z modular object
      */
    modular& operator*=(const modular& z) { return (*this) = (*this) * z; }

    /** Inplace Division.
      * self = self / z.
      * @return reference to self.
      * @param  z modular object
      */
    modular& operator/=(const modular& z) { return (*this) = (*this) / z; }

    //@} Inplace Arithmetic Operators

  private:
 
    friend ostream& operator<<(ostream&, const modular&);
    friend istream& operator>>(istream&, modular&);
    friend ostream& writeObject(ostream&, const modular&);

    /** Residue class representative.
      * Between 0 and (modulus - 1).
      */
    long residue;

    /** Prime modulus.
      * Static member means can only implement one field at a time.
      * This is a long, and not an unsigned long, because taking the remainer
      * of a negative long with an unsigned long gives results that are 
      * not always correct.
      */
    static long modulus;

    /** Reciprocal.
      * @return reference to 1 / self
      */
    modular& recip(void) const;

  }; // class modular

  modular& modular::recip(void) const
  {
    // The extended Euclidean algorithm
    int x, y, q, tx, ty, temp;
    x = modulus; y = residue;
    tx = 0; ty = 1;
    
#ifdef TRACE
    clog << "x = " << x << ", y = " << y << ", q = " << q 
         << ", tx = " << tx << ", ty = " << ty << endl;
#endif // TRACE

    while(y != 0)
    {
      // always: gcd(modulus,residue) = gcd(x,y)
      //	 sx*modulus + tx*residue = x
      //	 sy*modulus + ty*residue = y
      q = x / y; // integer quotient
      temp = y;  y  = x  - q*y;  x  = temp;
      temp = ty; ty = tx - q*ty; tx = temp;
    
#ifdef TRACE
    clog << "x = " << x << ", y = " << y << ", q = " << q 
         << ", tx = " << tx << ", ty = " << ty << endl;
#endif // TRACE

    }
    // now x = gcd(modulus,residue)
    modular* x_ptr = new modular(tx);
    return(*x_ptr);
  } // modular::recip()

  ostream& operator<<(ostream& os, const modular& z)
  {
    os << z.residue << " mod " << modular::modulus;
    return os;
  } // operator<<

  istream& operator>>(istream& is, modular &z)
  {
    long x;
    is >> x;
    z.put_residue(x);
    return is;
  } // operator>>

  long modular::modulus;  // declare static member

  // Specialization of LinBox functions

  /** @name Specialization of LinBox functions for class modular
   * These specializations allow the \Ref{unparam_field} and
   * \Ref{unparam_randIter} template classes wrap the class 
   * modular as a LinBox field.
   */
  //@{

  /** Conversion of field element to an integer.
   * This function assumes the output field element x has already been 
   * constructed, but that it is not already initialized.
   * This is done by converting the residue (obtained by get_residue)
   * to an integer through the use of static cast.
   * This, of course, assumes such static casts are possible.
   * @return reference to integer.
   * @param x reference to integer to contain output (reference returned).
   * @param y constant reference to field element.
   */
  template<> integer& 
  unparam_field<modular>::convert(integer& x, const modular& y) const
  { return x = static_cast<integer>(y.get_residue()); }

  /** Cardinality.
   * Return integer representing cardinality of the field.
   * In this specialization, returns the prime modulus of the field.
   * The modulus is obtained by get_modulus and then an integer is constructed
   * from this long.
   * This, of course, assumes such a constructor exists
   * @return integer representing cardinality of the field
   */
  template<> const integer& unparam_field<modular>::cardinality(void) const
  { return *(new integer(modular::get_modulus())); }

  /** Characteristic.
   * Return integer representing characteristic of the field.
   * In this specialization, returns the prime modulus of the field.
   * The modulus is obtained by get_modulus and then an integer is constructed
   * from this long.
   * This, of course, assumes such a constructor exists
   * @return integer representing characteristic of the field.
   */
  template<> const integer& unparam_field<modular>::characteristic(void) const
  { return *(new integer(modular::get_modulus())); }

  /** Faxpy apply method.
   * z = a*x + y.
   * @return reference to element z
   * @param z reference to element z
   * @param x constant reference to element x
   * @param y constant reference to element y
   */
  template<> unparam_field<modular>::element&
  faxpy<unparam_field<modular> >
    ::apply(unparam_field<modular>::element& z,
	    const unparam_field<modular>::element& x,
	    const unparam_field<modular>::element& y) const
  {
    long zint, 
         aint(_a.get_residue()), 
	 xint(x.get_residue()), 
	 yint(y.get_residue());
	 
    zint = aint * xint + yint;
    z.put_residue(zint);
    return z;
  }

  /** Faxpy inplace apply method.
   * y = a*x + y or y += a*x.
   * @return reference to element y
   * @param y reference to element y
   * @param x constant reference to element x
   */
  template<> unparam_field<modular>::element&
  faxpy<unparam_field<modular> >
    ::applyin(unparam_field<modular>::element& y,
	      const unparam_field<modular>::element& x) const
  {
    long zint,
  	 aint(_a.get_residue()), 
	 xint(x.get_residue()), 
	 yint(y.get_residue());
    zint = aint * xint + yint;
    y.put_residue(zint);
    return y;
  }

  //@} Specialization of LinBox functions

} // namespace LinBox

#endif // _MODULAR_FIELD_
