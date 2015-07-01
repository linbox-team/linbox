/* File: src/library/objects/field/fuzzy.h
 * Author: William Turner for the LinBox group
 */

#ifndef _FUZZY_FIELD_
#define _FUZZY_FIELD_

#include <iostream>
#include <cmath>
#include "LinBox/integer.h"
#include "LinBox/unparam_field.h"
#include "LinBox/faxpy.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

  /** Unparameterized field of "fuzzy" doubles.
    * Used with \Ref{unparam_field} to create implement a \Ref{LinBox}
    * unparamterized field of "fuzzy" doubles.
    * Field has static (non-negative) member to contain "fuzz value" of field.
    * Doubles within this fuzz value are considered to be equal.
    */
  class fuzzy
  {
  public:

    /** @name Object Management */
    //@{

    /** Default constructor.
      * Sets fuzz value to zero.
      */
    fuzzy () : value(0) {}
 
    /** Constructor from integer.
      * Sets residue to representative of residue class of integer.
      * @param  start prime integer modulus
      */
    fuzzy (integer start) : value(start) {}
 
    /** Copy constructor.
      * Copies residue.
      * @param  f fuzzy double object.
      */
    fuzzy (const fuzzy& f) { value = f.value; }

    /** Assignment Operator.
      * "Copies" into *this, which is an existing object.
      * If not explicitly defined, each member is assigned.
      * @return reference to self
      * @param  f fuzzy double object
      */
    fuzzy& operator=(const fuzzy& f)
    {
      if (this != &f ) value = f.value; // guard against self assignment
      return *this;
    }

    /** Retrieve double value.
      * @return double representing value
      */
    double get_value() const { return value; }

    /** Change value.
      * Changes value of fuzzy double to that of double.
      * @param  start double to which to change value of fuzzy double
      */
    void put_value(const double start) { value = start; }

    /** Retreive fuzz value.
      * @return fuzzy value
      */
    static double get_fuzz() { return fuzz; }
 
    /** Change fuzz value.
      * Changes fuzz value of field to double.
      * Negative values will be set to zero.
      * Static fuzz value means only one field can be implemented at a time.
      * @param  start double fuzz value
      */
    static void put_fuzz(const double start) 
    { 
      fuzz = start; 
      if (fuzz < 0) fuzz = 0; // do not allow negative fuzz values
    }

    //@} Object Management

    /** @name Arithmetic Operators */
    //@{

    /** Equality.
      * @return boolean true if elements are equal, false if not.
      * @param  f fuzzy double object.
      */
    bool operator==(const fuzzy& f) const
    { 
      // Allow equality in testing so that a fuzz value of 0 is allowed, which
      // enforeces strict equality and this class acts like the C++ doubles
      return ( ( value - f.value <= fuzz ) && ( f.value - value <= fuzz ) ); 
    }

    /** Addition.
      * @return self + f.
      * @param  f fuzzy double object
      */
    fuzzy& operator+(const fuzzy& f) const
    {
      fuzzy* x_ptr = new fuzzy;
      x_ptr->value = value + f.value;
      return *x_ptr;
    }

    /** Subtraction.
      * @return self - f.
      * @param  f fuzzy double object
      */
    fuzzy& operator-(const fuzzy& f) const
    {
      fuzzy* x_ptr = new fuzzy;
      x_ptr->value = value - f.value;
      return *x_ptr;
    }

    /** Multiplication.
      * @return self * f.
      * @param  f fuzzy double object
      */
    fuzzy& operator*(const fuzzy& f) const
    {
      fuzzy* x_ptr = new fuzzy;
      x_ptr->value = value * f.value;
      return *x_ptr;
    }

    /** Division.
      * @return self / f.
      * @param  f fuzzy double object
      */
    fuzzy& operator/(const fuzzy& f) const
    {
      fuzzy* x_ptr = new fuzzy;
      x_ptr->value = value / f.value;
      return *x_ptr;
    }

    /** Negation.
      * @return - self
      */
    fuzzy& operator-(void) const
    {
      fuzzy* x_ptr = new fuzzy;
      x_ptr->value = - value;
      return *x_ptr;
    }

    //@} Arithmetic Operators

    /** @name Inplace Arithmetic Operators */
    //@{

    /** Inplace Addition.
      * self = self + f.
      * @return reference to self.
      * @param  f fuzzy double object
      */
    fuzzy& operator+=(const fuzzy& f)
      { return (*this) = (*this) + f; }

    /** Inplace Subtraction.
      * self = self - f.
      * @return reference to self.
      * @param  f fuzzy double object
      */
    fuzzy& operator-=(const fuzzy& f)
      { return (*this) = (*this) - f; }

    /** Inplace Multiplication.
      * self = self * f.
      * @return reference to self.
      * @param  f fuzzy double object
      */
    fuzzy& operator*=(const fuzzy& f)
      { return (*this) = (*this) * f; }

    /** Inplace Division.
      * self = self / f.
      * @return reference to self.
      * @param  f fuzzy double object
      */
    fuzzy& operator/=(const fuzzy& f)
      { return (*this) = (*this) / f; }

    //@} Inplace Arithmetic Operators

  private:

    friend ostream& operator<<(ostream&, const fuzzy&);
    friend istream& operator>>(istream&, fuzzy&);

    /** Fuzzy double value.
      */
    double value;

    /** Fuzz value.
      * Static member means can only implement one field at a time.
      */
    static double fuzz;

  }; // class fuzzy

  double fuzzy::fuzz = 0;

  ostream& operator<<(ostream& os, const fuzzy& f)
  {
    os << f.value;
    return os;
  }

  istream& operator>>(istream& is, fuzzy &f)
  {
    is >> f.value;
    return is;
  }

  /** @name Specialization of LinBox functions for class fuzzy
   * These specializations allow the \Ref{unparam_field} and
   * \Ref{unparam_randIter} template classes wrap the class 
   * fuzzy as a LinBox field.
   */
  //@{

  /** Conversion of field element to an integer.
   * This function assumes the output field element x has already been 
   * constructed, but that it is not already initialized.
   * This is done by converting the value (obtained by get_value)
   * to an integer through the use of floor and static cast.
   * This, of course, assumes such static casts are possible.
   * @return reference to integer.
   * @param x reference to integer to contain output (reference returned).
   * @param y constant reference to field element.
   */
  template<> integer& 
  unparam_field<fuzzy>::convert(integer& x, const fuzzy& y) const
  { return x = static_cast<integer>(floor(y.get_value())); }

  //@} Specializations

} // namespace LinBox

#endif // _FUZZY_FIELD_
