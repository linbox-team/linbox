/* File: src/library/objects/field/param_fuzzy.h
 * Author: William Turner for the LinBox group
 */

#ifndef _PARAM_FUZZY_
#define _PARAM_FUZZY_

#include <iostream>
#include "LinBox/integer.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 

  // Forward declarations
  class param_fuzzy_randIter;

  /** Abstract parameterized field of "fuzzy" doubles.
   * Field has (non-static, non-negative) member to contain "fuzz value" of 
   * field.  Doubles within this fuzz value are considered to be equal.
   */
  class param_fuzzy
  {
  public:

    /** Element type.
     * It must meet the common object interface of elements as given in the
     * the archetype Element_archetype.
     */
    typedef double element;

    /** Random iterator generator type.
     * It must meet the common object interface of random element generators
     * as given in the the archetype RandIter_archetype.
     */
    typedef param_fuzzy_randIter randIter;

    /** @name Object Management
     */
    //@{
 
    /** Default constructor.
     */
    param_fuzzy(void) {}

    /** Constructor from an integer.
     * Sets the fuzz value of the field throug the static member of the 
     * element type.
     * @param value constant reference to double fuzz value
     */
    param_fuzzy(const double& value) : _fuzz(value) {}

    /** Copy constructor.
     * Constructs param_fuzzy object by copying the field.
     * This is required to allow field objects to be passed by value
     * into functions.
     * @param  F param_fuzzy object.
     */
    param_fuzzy(const param_fuzzy& F) : _fuzz(F._fuzz) {}
 
    /** Assignment operator.
     * @return reference to param_fuzzy object for self
     * @param F constant reference to param_fuzzy object
     */
    param_fuzzy& operator= (const param_fuzzy& F)
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
    { return x = static_cast<double>(y); }
 
   /** Conversion of field base element to a template class T.
     * This function assumes the output field base element x has already been
     * constructed, but that it is not already initialized.
     * @return reference to template class T.
     * @param x template class T to contain output (reference returned).
     * @param y constant field base element.
     */
    integer& convert(integer& x, const element& y) const
    { return x = static_cast<integer>(y); }
 
    /** Assignment of one field base element to another.
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    element& assign(element& x, const element& y) const
    { return x = y; }

    /** Cardinality.
     * Return integer representing cardinality of the domain.
     * Returns a non-negative integer for all domains with finite
     * cardinality, and returns -1 to signify a domain of infinite
     * cardinality.
     * @return integer representing cardinality of the domain
     */
    const integer& cardinality(void) const { return *(new integer(-1)); }
 
    /** Characteristic.
     * Return integer representing characteristic of the domain.
     * Returns a positive integer to all domains with finite characteristic,
     * and returns 0 to signify a domain of infinite characteristic.
     * @return integer representing characteristic of the domain.
     */
    const integer& characteristic(void) const { return *(new integer(0)); }

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
    { return ( ( x - y <= _fuzz ) && ( y - x <= _fuzz ) ); }

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
    { return x = y + z; }
 
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
    { return x = y - z; }
 
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
    { return x = y * z; }
 
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
    { return x = y / z; }
 
    /** Additive Inverse (Negation).
     * x = - y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    element& neg(element& x, const element& y) const
    { return x = - y; }
 
    /** Multiplicative Inverse.
     * x = 1 / y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    element& inv(element& x, const element& y) const
    { return x = static_cast<double>(1) / y; }

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
    { return ( ( x <= _fuzz ) && ( - x <= _fuzz ) ); }

    /** One equality.
     * Test if field base element is equal to one.
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return boolean true if equals one, false if not.
     * @param  x field base element.
     */
    bool isOne(const element& x) const
    { return ( ( x - 1 <= _fuzz ) && ( 1 - x <= _fuzz ) ); }

    /** Inplace Addition.
     * x += y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    element& addin(element& x, const element& y) const { return x += y; } 
    /** Inplace Subtraction.
     * x -= y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    element& subin(element& x, const element& y) const { return x -= y; }
 
    /** Inplace Multiplication.
     * x *= y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    element& mulin(element& x, const element& y) const { return x *= y; }
 
    /** Inplace Division.
     * x /= y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    element& divin(element& x, const element& y) const { return x /= y; }
 
    /** Inplace Additive Inverse (Inplace Negation).
     * x = - x
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     */
    element& negin(element& x) const { return x = -x; }
 
    /** Inplace Multiplicative Inverse.
     * x = 1 / x
     * This function assumes the field base elementhas already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     */
    element& invin(element& x) const
    { return x = static_cast<double>(1) / x; }

    //@} Inplace Arithmetic Operations

    /** @name Input/Output Operations */
    //@{

    /** Print field.
     * @return output stream to which field is written.
     * @param  os  output stream to which field is written.
     */
    ostream& write(ostream& os) const 
    { return os << " fuzz value " << _fuzz; }
 
    /** Read field.
     * @return input stream from which field is read.
     * @param  is  input stream from which field is read.
     */
    istream& read(istream& is) { return is >> _fuzz; }

    /** Print field base element.
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return output stream to which field base element is written.
     * @param  os  output stream to which field base element is written.
     * @param  x   field base element.
     */
    ostream& write(ostream& os, const element& x) const
    { return os << x; }
 
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
      return is; 
    }

    //@}

  private:

    /// Private static double for fuzz value
    double _fuzz;

  }; // class param_fuzzy

} // namespace LinBox

#include "LinBox/param_fuzzy_randiter.h"

#endif // _PARAM_FUZZY_
